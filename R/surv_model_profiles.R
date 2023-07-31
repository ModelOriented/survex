#' Helper functions for `model_profile.R`
#'
#' @param x an object containing calculated ceteris_paribus profiles
#' @param ... other parameters, ignored
#' @param variable_type character, either `"numerical"` or `"categorical"`, the type of variable to be calculated, if left `NULL` (default), both are calculated
#' @param variables a character vector containing names of variables to be explained
#' @param center logical, if the profiles should be centered before aggregations
#'
#' @return A data.frame with calculated results.
#'
#' @keywords internal
surv_aggregate_profiles <- function(x,
                                    ...,
                                    variable_type = NULL,
                                    groups = NULL,
                                    variables = NULL,
                                    center = FALSE) {
    all_profiles <- x$result
    class(all_profiles) <- "data.frame"

    all_profiles$`_ids_` <- factor(all_profiles$`_ids_`)


    # variables to use
    all_variables <-
        na.omit(as.character(unique(all_profiles$`_vname_`)))
    if (!is.null(variables)) {
        all_variables_intersect <- intersect(all_variables, variables)
        if (length(all_variables_intersect) == 0)
            stop(paste0(
                "parameter variables do not overlap with ",
                paste(all_variables, collapse = ", ")
            ))
        all_variables <- all_variables_intersect
    }

    if (!is.null(variable_type) && variable_type == "numerical") {
        all_profiles <-
            all_profiles[all_profiles$`_vtype_` == "numerical",]
    }

    if (!is.null(variable_type) && variable_type == "categorical") {
        all_profiles <-
            all_profiles[all_profiles$`_vtype_` == "categorical",]
    }

    all_variables <-
        intersect(all_variables, unique(all_profiles$`_vname_`))

    # select only suitable variables
    all_profiles <-
        all_profiles[all_profiles$`_vname_` %in% all_variables,]
    # create _x_
    tmp <- as.character(all_profiles$`_vname_`)
    for (viname in unique(tmp)) {
        all_profiles$`_x_`[tmp == viname] <-
            all_profiles[tmp == viname, viname]
    }

    if (!inherits(class(all_profiles), "data.frame")) {
        all_profiles <- as.data.frame(all_profiles)
    }

    # change x column to proper character values
    for (variable in all_variables) {
        if (variable %in% all_profiles[all_profiles$`_vtype_` == "categorical", "_vname_"])
            all_profiles[all_profiles$`_vname_` == variable,]$`_x_` <-
                as.character(apply(all_profiles[all_profiles$`_vname_` == variable,], 1, function(all_profiles)
                    all_profiles[all_profiles["_vname_"]]))
    }

    aggregated_profiles <-
        surv_aggregate_profiles_partial(all_profiles)
    class(aggregated_profiles) <-
        c(
            "aggregated_survival_profiles_explainer",
            "partial_dependence_survival_explainer",
            "data.frame"
        )

    return(aggregated_profiles)
}


surv_aggregate_profiles_partial <- function(all_profiles) {
    tmp <-
        all_profiles[, c("_vname_", "_vtype_", "_label_", "_x_", "_yhat_", "_times_")]
    aggregated_profiles <-
        aggregate(
            tmp$`_yhat_`,
            by = list(
                tmp$`_vname_`,
                tmp$`_vtype_`,
                tmp$`_label_`,
                tmp$`_x_`,
                tmp$`_times_`
            ),
            FUN = mean,
            na.rm = TRUE
        )
    colnames(aggregated_profiles) <-
        c("_vname_",  "_vtype_", "_label_", "_x_", "_times_", "_yhat_")
    aggregated_profiles$`_ids_` <- 0

    # for factors, keep proper order
    # as in https://github.com/ModelOriented/ingredients/issues/82
    if (!is.numeric(all_profiles$`_x_`)) {
        aggregated_profiles$`_x_` <-
            factor(aggregated_profiles$`_x_`,
                   levels = unique(all_profiles$`_x_`))
        aggregated_profiles <-
            aggregated_profiles[order(aggregated_profiles$`_x_`),]
    }

    aggregated_profiles
}


#' @param x an explainer object - model preprocessed by the `explain()` function
#' @param ... other parameters, ignored
#' @param data data used to create explanations
#' @param variables a character vector containing names of variables to be explained
#' @param categorical_variables character, a vector of names of additional variables which should be treated as categorical (factors are automatically treated as categorical variables). If it contains variable names not present in the `variables` argument, they will be added at the end.
#' @param grid_points maximum number of points for profile calculations. Note that the final number of points may be lower than grid_points. Will be passed to internal function. By default `51`.
#' @param center logical, if the profiles should be centered before aggregations
#'
#' @return A data.frame with calculated results.
#'
#' @keywords internal
surv_ale <- function(x,
                     ...,
                     data,
                     variables,
                     categorical_variables,
                     grid_points,
                     center = FALSE) {
    test_explainer(
        x,
        has_data = TRUE,
        has_survival = TRUE,
        function_name = "surv_ale"
    )

    if (is.null(variables))
        variables <- colnames(data)

    # change categorical_features to column names
    if (is.numeric(categorical_variables))
        categorical_variables <- colnames(data)[categorical_variables]
    additional_categorical_variables <- categorical_variables
    factor_variables <- colnames(data)[sapply(data, is.factor)]
    categorical_variables <-
        unique(c(additional_categorical_variables, factor_variables))

    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

    # Make predictions for original levels
    predictions_original <- predict_survival_function(model = model,
                                                      newdata = data,
                                                      times = times)
    mean_pred <- colMeans(predictions_original)

    profiles <- lapply(variables, function(variable) {
        X_lower <- X_upper <- data
        variable_values <- data[, variable]

        if (variable %in% categorical_variables) {
            if (!is.factor(variable_values)){
                data[, variable] <- as.factor(data[, variable])
                variable_values <- as.factor(variable_values)
            }
            levels_original <- levels(droplevels(variable_values))
            levels_n <- nlevels(droplevels(variable_values))
            if (inherits(variable_values, "ordered")) {
                level_order <- 1:levels_n
            } else {
                level_order <- order_levels(data, variable)
            }

            # The new order of the levels
            levels_ordered <- levels_original[level_order]

            # The feature with the levels in the new order
            x_ordered <-
                order(level_order)[as.numeric(droplevels(variable_values))]

            # Filter rows which are not already at maximum or minimum level values
            row_ind_increase <- (1:nrow(data))[x_ordered < levels_n]
            row_ind_decrease <- (1:nrow(data))[x_ordered > 1]
            X_lower[row_ind_decrease, variable] <-
                levels_ordered[x_ordered[row_ind_decrease] - 1]
            X_upper[row_ind_increase, variable] <-
                levels_ordered[x_ordered[row_ind_increase] + 1]

            # Make predictions for decreased levels (excluding minimum levels)
            predictions_lower <-
                predict_survival_function(model = model,
                                          newdata = X_lower[row_ind_decrease,],
                                          times = times)

            # Make predictions for increased levels (excluding maximum levels)
            predictions_upper <-
                predict_survival_function(model = model,
                                          newdata = X_upper[row_ind_increase,],
                                          times = times)

            d_increase <-
                predictions_upper - predictions_original[row_ind_increase,]
            d_decrease <-
                predictions_original[row_ind_decrease,] - predictions_lower
            prediction_deltas <- rbind(d_increase, d_decrease)
            colnames(prediction_deltas) <- times


            deltas <- data.frame(
                interval = rep(c(x_ordered[row_ind_increase],
                                 x_ordered[row_ind_decrease] - 1),
                               each = length(times)),
                time = rep(times, times = nrow(prediction_deltas)),
                yhat = c(t(prediction_deltas))
            )

            deltas <-
                aggregate(yhat ~ interval + time,
                          data = deltas,
                          FUN = mean)
            deltas1 <- deltas[deltas$interval == 1,]
            deltas1$yhat <- 0
            deltas$interval <- deltas$interval + 1
            deltas <- rbind(deltas, deltas1)
            deltas <- deltas[order(deltas$time, deltas$interval),]
            rownames(deltas) <- NULL
            deltas$yhat_cumsum <-
                ave(deltas$yhat, deltas$time, FUN = cumsum)

            x_count <- as.numeric(table(variable_values))
            x_prob <- x_count / sum(x_count)

            ale_means <-
                aggregate(
                    yhat_cumsum ~ time,
                    data = deltas,
                    FUN = function(x) {
                        sum(x * x_prob[level_order])
                    }
                )
            colnames(ale_means)[2] <- "ale0"

            ale_values <- merge(deltas,
                                ale_means,
                                all.x = TRUE,
                                by = "time")

            ale_values$ale <-
                ale_values$yhat_cumsum - ale_values$ale0
            ale_values$level <- levels_ordered[ale_values$interval]

            ale_values <-
                ale_values[order(ale_values$interval, ale_values$time),]
            if (!center){
                ale_values$ale <- ale_values$ale + mean_pred
            }
            return(
                data.frame(
                    `_vname_` = variable,
                    `_vtype_` = "categorical",
                    `_label_` = label,
                    `_x_` = ale_values$level,
                    `_times_` = ale_values$time,
                    `_yhat_` = ale_values$ale,
                    `_ids_` = 0,
                    check.names = FALSE
                )
            )

        } else {
            # Number of quantile points for determined by grid length
            quantile_vals <- as.numeric(quantile(
                variable_values,
                seq(0.01, 1, length.out = grid_points),
                type = 1
            ))

            # Quantile points vector
            quantile_vec <- c(min(variable_values), quantile_vals)
            quantile_vec <- unique(quantile_vec)

            quantile_df <-
                data.frame(id = 1:length(quantile_vec),
                           value = quantile_vec)

            # Match feature instances to quantile intervals
            interval_index <-
                findInterval(variable_values, quantile_vec, left.open = TRUE)

            # Points in interval 0 should be in interval 1
            interval_index[interval_index == 0] <- 1

            # Prepare datasets with upper and lower interval limits replacing original feature values
            X_lower[, variable] <- quantile_vec[interval_index]
            X_upper[, variable] <- quantile_vec[interval_index + 1]
            # Get survival predictions for instances of upper and lower interval limits
            predictions_lower <-
                predict_survival_function(model = model,
                                          newdata = X_lower,
                                          times = times)
            predictions_upper <-
                predict_survival_function(model = model,
                                          newdata = X_upper,
                                          times = times)

            # First order finite differences
            prediction_deltas <-
                predictions_upper - predictions_lower
            # Rename columns to timepoints for which predictions were made
            colnames(prediction_deltas) <- times

            deltas <- data.frame(
                x = rep(X_lower[, variable], each = length(times)),
                interval = rep(interval_index, each = length(times)),
                time = rep(times, times = nrow(data)),
                yhat = c(t(prediction_deltas))
            )

            deltas <-
                aggregate(yhat ~ interval + time,
                          data = deltas,
                          FUN = mean)
            deltas$yhat_cumsum <-
                ave(deltas$yhat, deltas$time, FUN = cumsum)
            interval_n <- as.numeric(table(interval_index))
            n <- sum(interval_n)

            ale_means <-
                aggregate(
                    yhat_cumsum ~ time,
                    data = deltas,
                    FUN = function(x) {
                        sum(((c(
                            0, x[1:(length(x) - 1)]
                        ) + x) / 2) * interval_n / n)
                    }
                )
            colnames(ale_means)[2] <- "ale0"

            # Centering the ALEs to obtain final ALE values
            ale_values <- merge(deltas,
                                ale_means,
                                all.x = TRUE,
                                by = "time")

            ale_values$ale <-
                ale_values$yhat_cumsum - ale_values$ale0
            ale_values$interval <- ale_values$interval + 1
            ale_values1 <-
                ale_values[seq(1, nrow(ale_values), length(quantile_vec) - 1), ]
            ale_values1$interval <- 1
            ale_values <- rbind(ale_values, ale_values1)

            ale_values <- merge(ale_values,
                                quantile_df,
                                by.x = "interval",
                                by.y = "id")
            ale_values <-
                ale_values[order(ale_values$interval, ale_values$time),]
            ale_values$ale <- ale_values$ale + mean_pred
            if (!center){
                ale_values$ale <- ale_values$ale + mean_pred
            }

            return(
                data.frame(
                    `_vname_` = variable,
                    `_vtype_` = "numerical",
                    `_label_` = label,
                    `_x_` = ale_values$value,
                    `_times_` = ale_values$time,
                    `_yhat_` = ale_values$ale,
                    `_ids_` = 0,
                    check.names = FALSE
                )
            )
        }

    })

    profiles <- do.call(rbind, profiles)
    class(profiles) <- c(
        "aggregated_survival_profiles_explainer",
        "accumulated_local_effects_survival_explainer",
        "data.frame"
    )
    return(profiles)
}
