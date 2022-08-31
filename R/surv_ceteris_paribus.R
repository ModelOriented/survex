#' Helper functions for `predict_profile.R`
#' @rdname surv_ceteris_paribus
#' @keywords internal
surv_ceteris_paribus <- function(x, ...) UseMethod("surv_ceteris_paribus", x)

#' Helper functions for `predict_profile.R`
#'
#' @rdname surv_ceteris_paribus
#'
#' @param x a model to be explained, preprocessed by the `explain` function
#' @param new_observation a new observation for which predictions need to be explained
#' @param variables character, names of the variables to be included in the calculations
#' @param categorical_variables character vector, names of variables that should be treated as categories (factors are included by default)
#' @param variable_splits named list of splits for variables, in most cases created with internal functions. If NULL then it will be calculated based on validation data available in the explainer
#' @param grid_points maximum number of points for profile calculations. Note that the final number of points may be lower than grid_points. Will be passed to internal function. By default `101`.
#' @param variable_splits_type character, decides how variable grids should be calculated. Use `"quantiles"` for percentiles or `"uniform"` (default) to get uniform grid of points.
#' @param ... other parameters, currently ignored
#'
#' @return A data.frame containing the result of the calculation.
#'
#' @keywords internal
surv_ceteris_paribus.surv_explainer <- function(x,
                                                    new_observation,
                                                    variables = NULL,
                                                    categorical_variables = NULL,
                                                    variable_splits = NULL,
                                                    grid_points = 101,
                                                    variable_splits_type = "uniform",
                                                    ...) {

    test_explainer(x, has_data = TRUE, has_survival = TRUE, has_y = TRUE, function_name = "ceteris_paribus_survival")

    data <- x$data
    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

    surv_ceteris_paribus.default(x = model,
                                     data = data,
                                     predict_survival_function = predict_survival_function,
                                     new_observation = new_observation,
                                     variables = variables,
                                     categorical_variables = categorical_variables,
                                     variable_splits = variable_splits,
                                     grid_points = grid_points,
                                     variable_splits_type = variable_splits_type,
                                     variable_splits_with_obs = TRUE,
                                     label = label,
                                     times = times,
                                     ...)

}

surv_ceteris_paribus.default <- function(x,
                                    data,
                                    predict_survival_function = NULL,
                                    new_observation,
                                    variables = NULL,
                                    categorical_variables = NULL,
                                    variable_splits = NULL,
                                    grid_points = 101,
                                    variable_splits_type = "uniform",
                                    variable_splits_with_obs = TRUE,
                                    label = NULL,
                                    times = times,
                                    ...) {


    if (is.data.frame(data)) {
        common_variables <- intersect(colnames(new_observation), colnames(data))
        new_observation <- new_observation[, common_variables, drop = FALSE]
        data <- data[, common_variables, drop = FALSE]
    }

    # change categorical_features to column names
    if (is.numeric(categorical_variables)) categorical_variables <- colnames(data)[categorical_variables]
    additional_categorical_variables <- categorical_variables
    factor_variables <- colnames(data)[sapply(data, is.factor)]
    categorical_variables <- unique(c(additional_categorical_variables, factor_variables))



    # calculate splits
    if (is.null(variable_splits)) {
        if (is.null(data))
            stop("The ceteris_paribus() function requires explainers created with specified 'data'.")
        if (is.null(variables))
            variables <- colnames(data)
        variable_splits <- calculate_variable_split(data, variables = variables,
                                                    categorical_variables = categorical_variables,
                                                    grid_points = grid_points,
                                                    variable_splits_type = variable_splits_type,
                                                    new_observation = if (variable_splits_with_obs) new_observation else NA)

        }

    profiles <- calculate_variable_survival_profile(new_observation,
                                                    variable_splits,
                                                    x,
                                                    predict_survival_function,
                                                    times)

    profiles$`_vtype_` <- ifelse(profiles$`_vname_` %in% categorical_variables, "categorical", "numerical")



    col_yhat <- grep(colnames(profiles), pattern = "^_yhat_")

    attr(profiles, "times") <- times
    attr(profiles, "observations") <- new_observation

    ret <- list(eval_times = times,
         variable_values = new_observation,
         result = cbind(profiles, `_label_` = label))

    class(ret) <- c("surv_ceteris_paribus", "list")

    ret
}


calculate_variable_split <- function(data, variables = colnames(data), categorical_variables = NULL,  grid_points = 101, variable_splits_type = "quantiles", new_observation = NA) {
    UseMethod("calculate_variable_split", data)
}

#' @importFrom stats na.omit quantile
calculate_variable_split.default <- function(data, variables = colnames(data), categorical_variables = NULL, grid_points = 101, variable_splits_type = "quantiles", new_observation = NA) {
    variable_splits <- lapply(variables, function(var) {
        selected_column <- na.omit(data[, var])

        if (!(var %in% categorical_variables)) {
            probs <- seq(0, 1, length.out = grid_points)
            if (variable_splits_type == "quantiles") {
                selected_splits <- unique(quantile(selected_column, probs = probs))
            } else {
                selected_splits <- seq(min(selected_column, na.rm = TRUE), max(selected_column, na.rm = TRUE), length.out = grid_points)
            }
            if (!any(is.na(new_observation)))
                selected_splits <- sort(unique(c(selected_splits, na.omit(new_observation[, var]))))
        } else {
            if (any(is.na(new_observation))) {
                selected_splits <- sort(unique(selected_column))
            } else {
                selected_splits <- sort(unique(rbind(data[, var, drop = FALSE],
                                                     new_observation[, var, drop = FALSE])[, 1]))
            }
        }
        selected_splits
    })
    names(variable_splits) <- variables
    variable_splits
}



calculate_variable_survival_profile <- function(data, variable_splits, model, predict_survival_function = NULL, times = NULL, ...) {
    UseMethod("calculate_variable_survival_profile")
}

calculate_variable_survival_profile.default <- function(data, variable_splits, model, predict_survival_function = NULL, times = NULL, ...) {


    variables <- names(variable_splits)
    prog <- progressr::progressor(along = 1:(length(variables)))
    profiles <- lapply(variables, function(variable) {
        split_points <- variable_splits[[variable]]


        if (is.null(rownames(data))) {
            ids <- rep(1:nrow(data), each = length(split_points)) # it never goes here, because null rownames are automatically setted to 1:n
        } else {
            ids <- rep(rownames(data), each = length(split_points))
        }



        new_data <- data[rep(1:nrow(data), each = length(split_points)), , drop = FALSE]
        new_data[, variable] <- rep(split_points, nrow(data))

        yhat <- c(t(predict_survival_function(model, new_data, times)))

        new_data <- data.frame(new_data[rep(seq_len(nrow(new_data)), each = length(times)), ],
                               `_times_` = rep(times, times = nrow(new_data)),
                               `_yhat_` = yhat,
                               `_vname_` = variable,
                               `_ids_` = ids,
                               check.names = FALSE)
        prog()
        new_data
    })

    profile <- do.call(rbind, profiles)
    class(profile) <- c("individual_variable_profile", class(profile))
    profile

}
