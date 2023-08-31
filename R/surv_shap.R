#' Helper functions for `predict_parts.R`
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param new_observation new observations for which predictions need to be explained
#' @param output_type a character, either `"survival"` or `"chf"`. Determines which type of prediction should be used for explanations.
#' @param ... additional parameters, passed to internal functions
#' @param y_true a two element numeric vector or matrix of one row and two columns, the first element being the true observed time and the second the status of the observation, used for plotting
#' @param calculation_method a character, either `"kernelshap"` for use of `kernelshap` library (providing faster Kernel SHAP with refinements) or `"exact_kernel"` for exact Kernel SHAP estimation
#' @param aggregation_method a character, either `"integral"`, `"integral_absolute"`, `"mean_absolute"`, `"max_absolute"`, or `"sum_of_squares"`
#'
#' @return A list, containing the calculated SurvSHAP(t) results in the `result` field
#'
#' @section References:
#' - \[1\] Krzyziński, Mateusz, et al. ["SurvSHAP(t): Time-dependent explanations of machine learning survival models."](https://www.sciencedirect.com/science/article/pii/S0950705122013302) Knowledge-Based Systems 262 (2023): 110234
#'
#' @keywords internal
surv_shap <- function(explainer,
                      new_observation,
                      output_type,
                      ...,
                      y_true = NULL,
                      calculation_method = "kernelshap",
                      aggregation_method = "integral") {
    stopifnot(
        "`y_true` must be either a matrix with one per observation in `new_observation` or a vector of length == 2" = ifelse(
            !is.null(y_true),
            ifelse(
                is.matrix(y_true),
                nrow(new_observation) == nrow(y_true),
                is.null(dim(y_true)) && length(y_true) == 2L
            ),
            TRUE
        )
    )

    test_explainer(explainer, "surv_shap", has_data = TRUE, has_y = TRUE, has_survival = TRUE)

    # make this code also work for 1-row matrix
    col_index <- which(colnames(new_observation) %in% colnames(explainer$data))
    if (is.matrix(new_observation) && nrow(new_observation) == 1) {
        new_observation <- as.matrix(t(new_observation[, col_index]))
    } else {
        new_observation <- new_observation[, col_index]
    }

    if (ncol(explainer$data) != ncol(new_observation)) {
        stop("New observation and data have different number of columns (variables)")
    }
    if (!is.null(y_true)) {
        if (is.matrix(y_true)) {
            # above, we have already checked that nrows of observations are
            # identical to nrows of y_true; thus we do not need to index
            # the first row here
            y_true_ind <- y_true[, 2]
            y_true_time <- y_true[, 1]
        } else {
            y_true_ind <- y_true[2]
            y_true_time <- y_true[1]
        }
    }

    # hack to use rf-model death times as explainer death times, as
    # treeshap::ranger_surv_fun.unify extracts survival times directly
    # from the ranger object for calculating the predictions
    if (calculation_method == "treeshap") {
        if (inherits(explainer$model, "ranger")) {
            explainer$times <- explainer$model$unique.death.times
        } else {
            stop("Calculation method `treeshap` is currently only implemented for `ranger`.")
        }
    }

    res <- list()
    res$eval_times <- explainer$times
    # to display final object correctly, when is.matrix(new_observation) == TRUE
    res$variable_values <- as.data.frame(new_observation)
    res$result <- switch(calculation_method,
        "exact_kernel" = use_exact_shap(explainer, new_observation, output_type, ...),
        "kernelshap" = use_kernelshap(explainer, new_observation, output_type, ...),
        stop("Only `exact_kernel` and `kernelshap` calculation methods are implemented")
    )

    if (!is.null(y_true)) res$y_true <- c(y_true_time = y_true_time, y_true_ind = y_true_ind)

    res$aggregate <- lapply(res$result, aggregate_surv_shap, method = aggregation_method, times = res$eval_times)

    if (nrow(new_observation) > 1) {
        class(res) <- "aggregated_surv_shap"
        # res$aggregation_method <- aggregation_method
        res$n_observations <- nrow(new_observation)
    } else {
        class(res) <- "surv_shap"
        res$result <- res$result[[1]]
        res$aggregate <- res$aggregate[[1]]
    }

    return(res)
}

use_exact_shap <- function(explainer, new_observation, output_type, observation_aggregation_method, ...) {
    shap_values <- sapply(
        X = as.character(seq_len(nrow(new_observation))),
        FUN = function(i) {
            as.data.frame(shap_kernel(explainer, new_observation[as.integer(i), ], output_type, ...))
        },
        USE.NAMES = TRUE,
        simplify = FALSE
    )

    return(shap_values)
}


shap_kernel <- function(explainer, new_observation, output_type, ...) {
    timestamps <- explainer$times
    p <- ncol(explainer$data)


    target_sf <- predict(explainer, new_observation, times = timestamps, output_type = output_type)
    sfs <- predict(explainer, explainer$data, times = timestamps, output_type = output_type)
    baseline_sf <- apply(sfs, 2, mean)


    permutations <- expand.grid(rep(list(0:1), p))
    kernel_weights <- generate_shap_kernel_weights(permutations, p)

    shap_values <- calculate_shap_values(
        explainer,
        explainer$model,
        baseline_sf,
        as.data.frame(explainer$data),
        permutations, kernel_weights,
        as.data.frame(new_observation),
        timestamps
    )



    shap_values <- as.data.frame(shap_values, row.names = colnames(explainer$data))
    colnames(shap_values) <- paste("t=", timestamps, sep = "")

    return(t(shap_values))
}

generate_shap_kernel_weights <- function(permutations, p) {
    apply(permutations, 1, function(row) {
        row <- as.numeric(row)
        num_available_variables <- sum(row != 0)

        if (num_available_variables == 0 || num_available_variables == p) {
            1e12
        } else {
            (p - 1) / (choose(p, num_available_variables) * num_available_variables * (p - num_available_variables))
        }
    })
}


calculate_shap_values <- function(explainer, model, avg_survival_function, data, simplified_inputs, shap_kernel_weights, new_observation, timestamps) {
    w <- shap_kernel_weights
    X <- as.matrix(simplified_inputs)

    R <- solve(t(sqrt(w) * X) %*% (sqrt(w) * X)) %*% t(X * w)

    y <- make_prediction_for_simplified_input(explainer, model, data, simplified_inputs, new_observation, timestamps)

    y <- sweep(
        y,
        2,
        avg_survival_function
    )

    R %*% y
}

make_prediction_for_simplified_input <- function(explainer, model, data, simplified_inputs, new_observation, timestamps) {
    preds <- apply(simplified_inputs, 1, function(row) {
        row <- as.logical(row)

        X_tmp <- data
        X_tmp[, row] <- new_observation[, row]
        colnames(X_tmp) <- colnames(data)

        colMeans(explainer$predict_survival_function(model, X_tmp, timestamps))
    })

    return(t(preds))
}

aggregate_surv_shap <- function(survshap, times, method, ...) {
    switch(method,
        "sum_of_squares" = return(apply(survshap, 2, function(x) sum(x^2))),
        "mean_absolute" = return(apply(survshap, 2, function(x) mean(abs(x)))),
        "max_absolute" = return(apply(survshap, 2, function(x) max(abs(x)))),
        "integral" = return(apply(survshap, 2, function(x) calculate_integral(x, times, normalization = "t_max"))),
        "integral_absolute" = return(apply(survshap, 2, function(x) calculate_integral(abs(x), times, normalization = "t_max"))),
        stop("aggregation_method has to be one of 'integral', 'integral_absolute', 'mean_absolute', 'max_absolute', or 'sum_of_squares'")
    )
}


use_kernelshap <- function(explainer, new_observation, output_type, observation_aggregation_method,  ...) {
    predfun <- function(model, newdata) {

        if (output_type == "survival"){
            explainer$predict_survival_function(
                model,
                newdata,
                times = explainer$times
            )
        } else {
            explainer$predict_cumulative_hazard_function(
                model,
                newdata,
                times = explainer$times
            )
        }

    }

    shap_values <- sapply(
        X = as.character(seq_len(nrow(new_observation))),
        FUN = function(i) {
            tmp_res <- kernelshap::kernelshap(
                object = explainer$model,
                X = new_observation[as.integer(i), ],
                bg_X = explainer$data,
                pred_fun = predfun,
                verbose = FALSE
            )
            tmp_shap_values <- data.frame(t(sapply(tmp_res$S, cbind)))
            colnames(tmp_shap_values) <- colnames(tmp_res$X)
            rownames(tmp_shap_values) <- paste("t=", explainer$times, sep = "")
            tmp_shap_values
        },
        USE.NAMES = TRUE,
        simplify = FALSE
    )

    return(shap_values)
}


#' @keywords internal
aggregate_shap_multiple_observations <- function(shap_res_list, feature_names, aggregation_function) {
    if (length(shap_res_list) > 1) {
        shap_res_list <- lapply(shap_res_list, function(x) {
            x$rn <- rownames(x)
            x
        })

        full_survshap_results <- do.call("rbind", shap_res_list)
        rownames(full_survshap_results) <- NULL

        # compute arithmetic mean for each time-point and feature across
        # multiple observations

        tmp_res <- aggregate(full_survshap_results[, !colnames(full_survshap_results) %in% c("rn")],
            by = list(full_survshap_results$rn),
            FUN = aggregation_function
        )
        rownames(tmp_res) <- tmp_res$Group.1
        ordering <- order(as.numeric(substring(rownames(tmp_res), 3)))

        tmp_res <- tmp_res[ordering, !colnames(tmp_res) %in% c("rn", "Group.1")]
    } else {
        # no aggregation required
        tmp_res <- shap_res_list[[1]]
    }
    shap_values <- tmp_res
    # transform to data.frame to make everything compatible with
    # previous code
    return(shap_values)
}
