#' @return An object of class `model_profile_2d_survival`. It is a list with the element `result` containing the results of the calculation.
#'
#' @export
model_profile_2d <- function(explainer,
                             variables = NULL,
                             N = 100,
                             ...,
                             categorical_variables = NULL,
                             grid_points = 51,
                             variable_splits_type = "uniform",
                             type = "partial",
                             output_type = "survival")
    UseMethod("model_profile_2d", explainer)

#' @export
model_profile_2d.surv_explainer <- function(explainer,
                                         variables = NULL,
                                         N = 100,
                                         ...,
                                         categorical_variables = NULL,
                                         grid_points = 51,
                                         variable_splits_type = "uniform",
                                         type = "partial",
                                         output_type = "survival"
                                         ) {
    variables <- unique(variables, categorical_variables)
    if (output_type != "survival") {
        stop("Currently only `survival` output type is implemented")
    }
    test_explainer(explainer,
                   "model_profile",
                   has_data = TRUE,
                   has_survival = TRUE)

    data <- explainer$data
    if (!is.null(N) && N < nrow(data)) {
        ndata <- data[sample(1:nrow(data), N), ]
    } else {
        ndata <- data
    }

    if (type == "partial") {
        result <- surv_pdp_2d(
            explainer,
            data = ndata,
            variables = variables,
            categorical_variables = categorical_variables,
            grid_points = grid_points,
            variable_splits_type = variable_splits_type,
            ...
        )
    } else if (type == "accumulated") {
        result <- surv_ale_2d(
            explainer,
            data = ndata,
            variables = variables,
            categorical_variables = categorical_variables,
            grid_points = grid_points,
            ...
        )
    } else {
        stop("Currently only `partial` and `accumulated` types are implemented")
    }

    ret <- list(
        result = result,
        eval_times = unique(result$`_times_`),
        variables = variables,
        type = type
    )
    class(ret) <- c("model_profile_2d_survival", "list")
    return(ret)
}

#' @keywords internal
surv_pdp_2d <- function(x,
                        data,
                        variables,
                        categorical_variables,
                        grid_points,
                        variable_splits_type,
                        ...) {
    data <- x$data
    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

    # change categorical_features to column names
    if (is.numeric(categorical_variables)) categorical_variables <- colnames(data)[categorical_variables]
    additional_categorical_variables <- categorical_variables
    factor_variables <- colnames(data)[sapply(data, is.factor)]
    categorical_variables <- unique(c(additional_categorical_variables, factor_variables))

    if (is.null(variables) | !is.list(variables) | !all(sapply(variables, length) == 2))
        stop("'variables' must be specified as a list of pairs (two-element vectors)")

    unique_variables <- unlist(variables)
    variable_splits <- calculate_variable_split(data,
                                                variables = unique_variables,
                                                categorical_variables = categorical_variables,
                                                grid_points = grid_points,
                                                variable_splits_type = variable_splits_type)

    profiles <- lapply(variables, FUN = function(variables_pair){
        var1 <- variables_pair[1]
        var2 <- variables_pair[2]
        expanded_data <- merge(variable_splits[[var1]], data[,!colnames(data) %in% variables_pair])
        names(expanded_data)[colnames(expanded_data) == "x"] <- var1
        expanded_data <- merge(variable_splits[[var2]], expanded_data)
        names(expanded_data)[colnames(expanded_data) == "x"] <- var2
        expanded_data <- expanded_data[,colnames(data)]

        predictions <- predict_survival_function(model = model,
                                       newdata = expanded_data,
                                       times = times)
        res <- data.frame(
            "_v1name_" = var1,
            "_v2name_" = var2,
            "_v1type_" = ifelse(var1 %in% categorical_variables, "categorical", "numerical"),
            "_v2type_" = ifelse(var2 %in% categorical_variables, "categorical", "numerical"),
            "_v1value_" = rep(expanded_data[,var1], each=length(times)),
            "_v2value_" = rep(expanded_data[,var2], each=length(times)),
            "_times_" = rep(times, nrow(expanded_data)),
            "_yhat_" = c(t(predictions)),
            "_label_" = label,
            check.names = FALSE
        )
        return(aggregate(`_yhat_`~., data = res, FUN=mean))
        })

    profiles <- do.call(rbind, profiles)
    profiles
}

