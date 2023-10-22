#' Helper functions for `predict_parts.R`
#'
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param new_observation a new observation for which predictions need to be explained
#' @param ... additional parameters, passed to internal functions
#' @param N a positive integer, number of observations generated in the neighbourhood
#' @param distance_metric character, name of the distance metric to be used, only `"euclidean"` is implemented
#' @param kernel_width a numeric, parameter used for calculating weights, by default it's `sqrt(ncol(data)*0.75)`
#' @param sampling_method character, name of the method of generating neighbourhood, only `"gaussian"` is implemented
#' @param sample_around_instance logical, if the neighbourhood should be generated with the new observation as the center (default), or should the mean of the whole dataset be used as the center
#' @param max_iter a numeric, maximal number of iteration for the optimization problem
#' @param categorical_variables character vector, names of variables that should be treated as categories (factors are included by default)
#' @param k a small positive number > 1, added to chf before taking log, so that weigths aren't negative
#'
#' @return A list, with the SurvLIME result in the `$result` field.
#'
#' @section References:
#' - \[1\] Kovalev, Maxim S., et al. ["SurvLIME: A method for explaining machine learning survival models."](https://www.sciencedirect.com/science/article/pii/S0950705120304044?casa_token=6e9cyk_ji3AAAAAA:tbqo33MsZvNC9nrSGabZdLfPtZTsvsvZTHYQCM2aEhumLI5D46U7ovhr37EaYUhmKZrw45JzDhg) Knowledge-Based Systems 203 (2020): 106164.
#'
#' @keywords internal
#' @importFrom stats optim
surv_lime <- function(explainer, new_observation,
                      ...,
                      N = 100,
                      distance_metric = "euclidean",
                      kernel_width = NULL,
                      sampling_method = "gaussian",
                      sample_around_instance = TRUE,
                      max_iter = 10000,
                      categorical_variables = NULL,
                      k = 1 + 1e-4) {
    test_explainer(explainer, "surv_lime", has_data = TRUE, has_y = TRUE, has_chf = TRUE)
    new_observation <- new_observation[, colnames(new_observation) %in% colnames(explainer$data)]
    if (ncol(explainer$data) != ncol(new_observation)) stop("New observation and data have different number of columns (variables)")

    predicted_sf <- explainer$predict_survival_function(explainer$model, new_observation, explainer$times)

    neighbourhood <- generate_neighbourhood(
        explainer$data,
        new_observation,
        N,
        categorical_variables,
        sampling_method,
        sample_around_instance
    )


    sc <- attr(neighbourhood$inverse, "scaled:scale")
    me <- attr(neighbourhood$inverse, "scaled:center")

    scaled_data <- neighbourhood$data


    scaled_data[, names(sc)] <- (scaled_data[, names(sc)] - me[col(scaled_data[, names(sc)])]) / sc[col(scaled_data[, names(sc)])]
    dist <- function(x, y) sqrt(sum((x - y)^2))

    distances <- apply(scaled_data, 1, dist, scaled_data[1, ])

    if (is.null(kernel_width)) kernel_width <- sqrt(ncol(scaled_data) * 0.75)

    weights <- sqrt(exp(-(distances^2) / (kernel_width^2)))
    na_est <- survival::basehaz(survival::coxph(explainer$y ~ 1))

    model_chfs <- explainer$predict_cumulative_hazard_function(explainer$model, neighbourhood$inverse, na_est$time) + k
    log_chfs <- log(model_chfs)
    weights_v <- model_chfs / log_chfs
    t_diffs <- c(diff(na_est$time), 1e-32)


    loss <- function(beta) {
        multiplied <- as.matrix(neighbourhood$inverse_ohe) %*% as.matrix(beta)
        multiplied <- multiplied[, rep(1, times = ncol(model_chfs))]

        log_na_est <- t(matrix(rep(log(na_est$hazard + k), N), ncol = N))
        sum(
            weights *
                rowSums(
                    (weights_v^2) *
                        ((log_chfs - log_na_est - multiplied)^2)
                        * t_diffs
                )
        )
    }

    var_values <- neighbourhood$inverse_ohe[1, ]

    res <- optim(rep(0, times = ncol(neighbourhood$inverse_ohe)), loss)

    beta <- data.frame(t(res$par))
    names(beta) <- colnames(neighbourhood$inverse_ohe)
    ret_list <- list(
        result = beta,
        variable_values = var_values,
        black_box_sf_times = explainer$times,
        black_box_sf = as.numeric(explainer$predict_survival_function(explainer$model, new_observation, explainer$times)),
        expl_sf_times = na_est$time,
        expl_sf = exp(-na_est$hazard * exp(sum(var_values * res$par)))
    )

    class(ret_list) <- c("surv_lime", class(ret_list))
    return(ret_list)
}


#' @importFrom stats rnorm model.matrix as.formula
generate_neighbourhood <- function(data_org,
                                   data_row,
                                   n_samples = 100,
                                   categorical_variables = NULL,
                                   sampling_method = "gaussian",
                                   sample_around_instance = TRUE) {

    # change categorical_variables to column names
    if (is.numeric(categorical_variables)) categorical_variables <- colnames(data_org)[categorical_variables]
    additional_categorical_variables <- categorical_variables
    factor_variables <- colnames(data_org)[sapply(data_org, is.factor)]
    categorical_variables <- unique(c(additional_categorical_variables, factor_variables))
    data_row <- data_row[colnames(data_org)]

    feature_frequencies <- list(length(categorical_variables))
    scaled_data <- scale(data_org[, !colnames(data_org) %in% categorical_variables])

    for (feature in categorical_variables) {
        column <- data_org[, feature]

        if (!is.factor(column)) column <- as.factor(column)
        feature_count <- summary(column)
        frequencies <- feature_count / sum(feature_count)
        feature_frequencies[[feature]] <- frequencies
    }

    n_col <- ncol(data_org[, !colnames(data_org) %in% categorical_variables])
    sc <- attr(scaled_data, "scaled:scale")
    me <- attr(scaled_data, "scaled:center")

    data <- switch(sampling_method,
        "gaussian" = matrix(rnorm(n_samples * n_col), nrow = n_samples, ncol = n_col),
        stop("Only `gaussian` sampling_method is implemented")
    )

    if (sample_around_instance) {
        to_add <- data_row[, !colnames(data_row) %in% categorical_variables]
        data <- data %*% diag(sc) + to_add[col(data)]
    } else {
        data <- data %*% diag(sc) + me[col(data)]
    }

    data <- as.data.frame(data)
    colnames(data) <- names(me)
    inverse <- data

    for (feature in categorical_variables) {
        values <- names(feature_frequencies[[feature]])
        freqs <- feature_frequencies[[feature]]
        inverse_column <- sample(values, n_samples, replace = TRUE, prob = freqs)
        if (feature %in% additional_categorical_variables) {
            inverse[, feature] <- as.numeric(inverse_column)
            data_row[, feature] <- as.numeric(data_row[, feature])
        } else {
            inverse[, feature] <- inverse_column
            data_row[, feature] <- as.character(data_row[, feature])
        }
        data[, feature] <- as.numeric(inverse_column == rep(data_row[feature], n_samples))
        data[1, feature] <- 1
    }


    inverse[1, ] <- data_row[1, colnames(inverse)]
    inverse <- inverse[, colnames(data_row)]
    attr(inverse, "scaled:scale") <- attr(scaled_data, "scaled:scale")
    attr(inverse, "scaled:center") <- attr(scaled_data, "scaled:center")

    data <- data[, colnames(data_row)]

    if (length(categorical_variables) > 0) {
        inverse_as_factor <- inverse
        inverse_as_factor[additional_categorical_variables] <-
            lapply(inverse_as_factor[additional_categorical_variables], as.factor)
        expr <- paste0("~", paste(categorical_variables, collapse = "+"))
        categorical_matrix <- model.matrix(as.formula(expr), data = inverse_as_factor)[, -1]
        inverse_ohe <- cbind(inverse, categorical_matrix)
        inverse_ohe[, categorical_variables] <- NULL
    } else {
        inverse_ohe <- inverse
    }

    list(data = data, inverse = inverse, inverse_ohe = inverse_ohe)
}
