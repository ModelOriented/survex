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

    if (is.null(variables) | !is.list(variables) | !all(sapply(variables, length) == 2))
        stop("'variables' must be specified as a list of pairs (two-element vectors)")

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

    # change categorical_features to column names
    if (is.numeric(categorical_variables)) categorical_variables <- colnames(data)[categorical_variables]
    additional_categorical_variables <- categorical_variables
    factor_variables <- colnames(data)[sapply(data, is.factor)]
    categorical_variables <- unique(c(additional_categorical_variables, factor_variables))

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

surv_pdp_2d <- function(x,
                        data,
                        variables,
                        categorical_variables,
                        grid_points,
                        variable_splits_type,
                        ...) {
    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

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

surv_ale_2d <- function(x,
                        data,
                        variables,
                        categorical_variables,
                        grid_points,
                        ...){
    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

    predictions_original <- predict_survival_function(model = model,
                                                      newdata = data,
                                                      times = times)
    mean_pred <- colMeans(predictions_original)


    profiles <- lapply(variables, FUN = function(variables_pair){
        var1 <- variables_pair[1]
        var2 <- variables_pair[2]

        if (all(!variables_pair %in% categorical_variables)){
            surv_ale_2d_num_num(
                model,
                data,
                predict_survival_function,
                times,
                grid_points,
                var1,
                var2,
            )
        }

    })


}


surv_ale_2d_num_num <- function(model,
                                data,
                                predict_survival_function,
                                times,
                                grid_points,
                                var1,
                                var2){

    # Number of quantile points for determined by grid length
    quantile_vals1 <- as.numeric(quantile(data[, var1],
                                          seq(0.01, 1, length.out = grid_points),
                                          type = 1))
    quantile_vals2 <- as.numeric(quantile(data[, var2],
                                          seq(0.01, 1, length.out = grid_points),
                                          type = 1))

    quantile_vec1 <- unique(c(min(data[, var1]), quantile_vals1))
    quantile_vec2 <- unique(c(min(data[, var2]), quantile_vals2))

    data <- data[(data[, var1] <= max(quantile_vec1)) &
                           (data[, var1] >= min(quantile_vec1)) &
                           (data[, var2] <= max(quantile_vec2)) &
                           (data[, var2] >= min(quantile_vec2)), ]

    # Matching instances to the grids of both features
    interval_index1 <- findInterval(data[, var1], quantile_vec1, left.open = TRUE)
    interval_index2 <- findInterval(data[, var2], quantile_vec2, left.open = TRUE)

    interval_index1[interval_index1 == 0] <- 1
    interval_index2[interval_index2 == 0] <- 1

    X_low1_low2 <- X_up1_low2 <- X_low1_up2 <- X_up1_up2 <- data
    X_low1_low2[, c(var1, var2)] <- cbind(quantile_vec1[interval_index1],
                                          quantile_vec2[interval_index2])
    X_up1_low2[, c(var1, var2)] <- cbind(quantile_vec1[interval_index1 + 1],
                                         quantile_vec2[interval_index2])
    X_low1_up2[, c(var1, var2)] <- cbind(quantile_vec1[interval_index1],
                                         quantile_vec2[interval_index2 + 1])
    X_up1_up2[, c(var1, var2)] <- cbind(quantile_vec1[interval_index1 + 1],
                                        quantile_vec2[interval_index2 + 1])

    y_hat_11 <- predict_survival_function(model = model, newdata = X_low1_low2, times = times)
    y_hat_21 <- predict_survival_function(model = model, newdata = X_up1_low2, times = times)
    y_hat_12 <- predict_survival_function(model = model, newdata = X_low1_up2, times = times)
    y_hat_22 <- predict_survival_function(model = model, newdata = X_up1_up2, times = times)

    prediction_deltas <- (y_hat_22 - y_hat_21) - (y_hat_12 - y_hat_11)
    colnames(prediction_deltas) <- times

    deltas <- data.frame(
        interval1 = rep(interval_index1, each = length(times)),
        interval2 = rep(interval_index2, each = length(times)),
        time = rep(times, times = nrow(data)),
        yhat = c(t(prediction_deltas))
    )

    deltas <- aggregate(`yhat`~., data = deltas, FUN=mean)
    interval_grid <- expand.grid(
        interval1 = unique(deltas$interval1),
        interval2 = unique(deltas$interval2),
        time = times
    )
    deltas <- merge(deltas,
                    interval_grid,
                    on = c("interval1", "interval2"),
                    all.y = TRUE)

    deltas$yhat_cumsum <- ave(deltas$yhat, deltas$time, deltas$interval1, FUN = function(x) cumsum(ifelse(is.na(x), 0, x)))
    deltas$yhat_cumsum <- ave(deltas$yhat_cumsum, deltas$time, deltas$interval2, FUN = function(x) cumsum(ifelse(is.na(x), 0, x)))

    cell_counts <- as.matrix(table(interval_index1, interval_index2))
    cell_counts_df <- as.data.frame(as.table(cell_counts))
    colnames(cell_counts_df) <- c("interval1", "interval2", "count")
    cell_counts_df$interval1 <- as.numeric(as.character(cell_counts_df$interval1))
    cell_counts_df$interval2 <- as.numeric(as.character(cell_counts_df$interval2))
    ale <- merge(deltas, cell_counts_df, on = c("interval1", "interval2"), all.x = TRUE)


    # Computing the first-order effect of feature 1
    ale1 <- data.frame()  # Initialize an empty data frame

    res <- copy(ale)
    res$yhat_diff <- ave(ale$yhat_cumsum,
            list(ale$interval2, ale$time),
                  FUN = function(x) c(x[1], x[-1] - x[-length(x)]))

    res$ale1

    res$ale1 <- ave(res$yhat_diff,
        list(res$interval1, res$time),
        FUN = function(x) c((x[-length(x)] + x[-1]) / 2, 0))

    # sub_res <- res[, c("time", "interval1", "count", "yhat_diff")]
    #
    # # Step 2: Calculate the numerator and denominator for each group using ave
    # sub_res$numerator <- with(sub_res, ave(count[-1] * (yhat_diff[-nrow(sub_res)] + yhat_diff[-1]) / 2,
    #                                        time, interval1, FUN = sum))
    # sub_res$denominator <- with(sub_res, ave(count[-1], time, interval1, FUN = sum))


}
