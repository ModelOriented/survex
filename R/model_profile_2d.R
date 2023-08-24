#' Dataset Level 2-Dimensional Variable Profile for Survival Models
#'
#' This function calculates explanations on a dataset level that help explore model response as a function of selected pairs of variables.
#' The explanations are calculated as an extension of Partial Dependence Profiles or Accumulated Local Effects with the inclusion of the time dimension.
#'
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param variables list of character vectors of length 2, names of pairs of variables to be explained
#' @param N number of observations used for the calculation of aggregated profiles. By default `100`. If `NULL` all observations are used.
#' @param categorical_variables character, a vector of names of additional variables which should be treated as categorical (factors are automatically treated as categorical variables). If it contains variable names not present in the `variables` argument, they will be added at the end.
#' @param grid_points maximum number of points for profile calculations. Note that the final number of points may be lower than grid_points. Will be passed to internal function. By default `25`.
#' @param center logical, should profiles be centered around the average prediction
#' @param variable_splits_type character, decides how variable grids should be calculated. Use `"quantiles"` for quantiles or `"uniform"` (default) to get uniform grid of points. Used only if `type = "partial"`.
#' @param type the type of variable profile, `"partial"` for Partial Dependence or `"accumulated"` for Accumulated Local Effects
#' @param output_type either `"survival"` or `"risk"` the type of survival model output that should be considered for explanations. Currently only `"survival"` is available.
#'
#' @return An object of class `model_profile_2d_survival`. It is a list with the element `result` containing the results of the calculation.
#'
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' cph_model_profile_2d <- model_profile_2d(cph_exp,
#'     variables = list(c("age", "celltype"))
#' )
#' head(cph_model_profile_2d$result)
#' plot(cph_model_profile_2d)
#'
#' cph_model_profile_2d_ale <- model_profile_2d(cph_exp,
#'     variables = list(c("age", "karno")),
#'     type = "accumulated"
#' )
#' head(cph_model_profile_2d_ale$result)
#' plot(cph_model_profile_2d_ale)
#' }
#'
#' @rdname model_profile_2d.surv_explainer
#' @export
model_profile_2d <- function(explainer,
                             variables = NULL,
                             N = 100,
                             categorical_variables = NULL,
                             grid_points = 25,
                             center = FALSE,
                             variable_splits_type = "uniform",
                             type = "partial",
                             output_type = "survival") {
    UseMethod("model_profile_2d", explainer)
}


#' @rdname model_profile_2d.surv_explainer
#' @export
model_profile_2d.surv_explainer <- function(explainer,
                                            variables = NULL,
                                            N = 100,
                                            categorical_variables = NULL,
                                            grid_points = 25,
                                            center = FALSE,
                                            variable_splits_type = "uniform",
                                            type = "partial",
                                            output_type = "survival") {
    if (is.null(variables) || !is.list(variables) || !all(sapply(variables, length) == 2)) {
        stop("'variables' must be specified as a list of pairs (two-element vectors)")
    }

    if (output_type != "survival") {
        stop("Currently only `survival` output type is implemented")
    }
    test_explainer(explainer, "model_profile", has_data = TRUE, has_survival = TRUE)

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
            center = center
        )
    } else if (type == "accumulated") {
        result <- surv_ale_2d(
            explainer,
            data = ndata,
            variables = variables,
            categorical_variables = categorical_variables,
            grid_points = grid_points,
            center = center
        )
    } else {
        stop("Currently only `partial` and `accumulated` types are implemented")
    }

    ret <- list(
        result = result,
        eval_times = unique(result$`_times_`),
        variables = variables,
        type = type,
        median_survival_time = explainer$median_survival_time
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
                        center) {
    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

    unique_variables <- unlist(variables)
    variable_splits <- calculate_variable_split(data,
        variables = unique_variables,
        categorical_variables = categorical_variables,
        grid_points = grid_points,
        variable_splits_type = variable_splits_type
    )

    profiles <- lapply(variables, FUN = function(variables_pair) {
        var1 <- variables_pair[1]
        var2 <- variables_pair[2]
        expanded_data <- merge(variable_splits[[var1]], data[, !colnames(data) %in% variables_pair])
        names(expanded_data)[colnames(expanded_data) == "x"] <- var1
        expanded_data <- merge(variable_splits[[var2]], expanded_data)
        names(expanded_data)[colnames(expanded_data) == "x"] <- var2
        expanded_data <- expanded_data[, colnames(data)]


        predictions_original <- predict_survival_function(
            model = model,
            newdata = data,
            times = times
        )
        mean_pred <- colMeans(predictions_original)

        predictions <- predict_survival_function(
            model = model,
            newdata = expanded_data,
            times = times
        )

        preds <- c(t(predictions))
        if (center) {
            preds <- preds - mean_pred
        }

        res <- data.frame(
            "_v1name_" = var1,
            "_v2name_" = var2,
            "_v1type_" = ifelse(var1 %in% categorical_variables, "categorical", "numerical"),
            "_v2type_" = ifelse(var2 %in% categorical_variables, "categorical", "numerical"),
            "_v1value_" = as.character(rep(expanded_data[, var1], each = length(times))),
            "_v2value_" = as.character(rep(expanded_data[, var2], each = length(times))),
            "_times_" = rep(times, nrow(expanded_data)),
            "_yhat_" = preds,
            "_label_" = label,
            check.names = FALSE
        )
        return(aggregate(`_yhat_` ~ ., data = res, FUN = mean))
    })

    profiles <- do.call(rbind, profiles)
    profiles
}

surv_ale_2d <- function(x,
                        data,
                        variables,
                        categorical_variables,
                        grid_points,
                        center) {
    model <- x$model
    label <- x$label
    predict_survival_function <- x$predict_survival_function
    times <- x$times

    predictions_original <- predict_survival_function(
        model = model,
        newdata = data,
        times = times
    )
    mean_pred <- colMeans(predictions_original)


    profiles <- lapply(variables, FUN = function(variables_pair) {
        var1 <- variables_pair[1]
        var2 <- variables_pair[2]

        if (all(!variables_pair %in% categorical_variables)) {
            surv_ale_2d_num_num(
                model,
                data,
                label,
                predict_survival_function,
                times,
                grid_points,
                var1,
                var2,
                mean_pred,
                center
            )
        } else {
            stop("Currently 2D ALE are implemented only for pairs of numerical variables")
        }
    })

    profiles <- do.call(rbind, profiles)
    profiles
}

#' @importFrom stats reshape
surv_ale_2d_num_num <- function(model,
                                data,
                                label,
                                predict_survival_function,
                                times,
                                grid_points,
                                var1,
                                var2,
                                mean_pred,
                                center) {

    # Number of quantile points for determined by grid length
    quantile_vals1 <- as.numeric(quantile(data[, var1],
        seq(0.01, 1, length.out = grid_points),
        type = 1
    ))
    quantile_vals2 <- as.numeric(quantile(data[, var2],
        seq(0.01, 1, length.out = grid_points),
        type = 1
    ))

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
    X_low1_low2[, c(var1, var2)] <- cbind(
        quantile_vec1[interval_index1],
        quantile_vec2[interval_index2]
    )
    X_up1_low2[, c(var1, var2)] <- cbind(
        quantile_vec1[interval_index1 + 1],
        quantile_vec2[interval_index2]
    )
    X_low1_up2[, c(var1, var2)] <- cbind(
        quantile_vec1[interval_index1],
        quantile_vec2[interval_index2 + 1]
    )
    X_up1_up2[, c(var1, var2)] <- cbind(
        quantile_vec1[interval_index1 + 1],
        quantile_vec2[interval_index2 + 1]
    )

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

    deltas <- aggregate(`yhat` ~ ., data = deltas, FUN = mean)
    interval_grid <- expand.grid(
        interval1 = c(0, sort(unique(deltas$interval1))),
        interval2 = c(0, sort(unique(deltas$interval2))),
        time = times
    )
    deltas <- merge(deltas,
        interval_grid,
        on = c("interval1", "interval2"),
        all.y = TRUE
    )

    deltas$yhat_cumsum <- ave(deltas$yhat, deltas$time, deltas$interval1, FUN = function(x) cumsum(ifelse(is.na(x), 0, x)))
    deltas$yhat_cumsum <- ave(deltas$yhat_cumsum, deltas$time, deltas$interval2, FUN = function(x) cumsum(ifelse(is.na(x), 0, x)))

    interval_index1 <- factor(interval_index1, sort(unique(interval_index1)))
    interval_index2 <- factor(interval_index2, sort(unique(interval_index2)))
    cell_counts <- as.matrix(table(interval_index1, interval_index2))
    cell_counts_df <- as.data.frame(as.table(cell_counts))
    colnames(cell_counts_df) <- c("interval1", "interval2", "count")
    cell_counts_df$interval1 <- as.numeric(as.character(cell_counts_df$interval1))
    cell_counts_df$interval2 <- as.numeric(as.character(cell_counts_df$interval2))
    ale <- merge(deltas, cell_counts_df, on = c("interval1", "interval2"), all.x = TRUE)
    ale <- ale[order(ale$interval1, ale$interval2), ]

    # Computing the first-order effect of feature 1
    res <- ale
    res$yhat_diff <- ave(ale$yhat_cumsum,
        list(ale$interval2, ale$time),
        FUN = function(x) c(x[1], diff(x))
    )

    ale1 <- do.call("rbind", lapply(sort(unique(res$interval1)), function(x) {
        counts <- res[res$interval1 == x & res$time == times[1], "count"]
        aggregate(yhat_diff ~ time,
            data = res[res$interval1 == x, ],
            FUN = function(vals) {
                sum(counts[-1] * (vals[-length(vals)] + vals[-1]) / 2 / sum(counts[-1]))
            }
        )
    }))

    ale1$interval1 <- rep(sort(unique(res$interval1)), each = length(times))
    ale1$yhat_diff[is.na(ale1$yhat_diff)] <- 0
    ale1$ale1 <- ave(ale1$yhat_diff, ale1$time, FUN = cumsum)


    # Computing the first-order effect of feature 2
    res <- ale
    res$yhat_diff <- ave(ale$yhat_cumsum,
        list(ale$interval1, ale$time),
        FUN = function(x) c(x[1], diff(x))
    )

    ale2 <- do.call("rbind", lapply(sort(unique(res$interval2)), function(x) {
        counts <- res[res$interval2 == x & res$time == times[1], "count"]
        aggregate(yhat_diff ~ time,
            data = res[res$interval2 == x, ],
            FUN = function(vals) {
                sum(counts[-1] * (vals[-length(vals)] + vals[-1]) / 2 / sum(counts[-1]))
            }
        )
    }))

    ale2$interval2 <- rep(sort(unique(res$interval2)), each = length(times))
    ale2$yhat_diff[is.na(ale2$yhat_diff)] <- 0
    ale2$ale2 <- ave(ale2$yhat_diff, ale2$time, FUN = cumsum)


    fJ0 <- unlist(lapply(times, function(time) {
        ale_time <- ale[ale$time == time, ]
        ale1_time <- ale1[ale1$time == time, ]
        ale2_time <- ale2[ale2$time == time, ]

        ale_time <- ale_time[c("interval1", "interval2", "yhat_cumsum")]
        dd <- reshape(ale_time,
            idvar = "interval1",
            timevar = "interval2",
            direction = "wide"
        )[, -1]
        rownames(dd) <- unique(ale_time$interval1)
        colnames(dd) <- unique(ale_time$interval2)

        dd <- dd - outer(ale1_time$ale1, rep(1, nrow(ale2_time))) -
            outer(rep(1, nrow(ale1_time)), ale2_time$ale2)
        sum(cell_counts * (dd[1:(nrow(dd) - 1), 1:(ncol(dd) - 1)] +
            dd[1:(nrow(dd) - 1), 2:ncol(dd)] +
            dd[2:nrow(dd), 1:(ncol(dd) - 1)] +
            dd[2:nrow(dd), 2:ncol(dd)]) / 4, na.rm = TRUE) / sum(cell_counts)
    }))

    fJ0 <- data.frame("fJ0" = fJ0, time = times)
    ale <- merge(ale, fJ0, by = c("time"))
    ale <- merge(ale, ale1, by = c("time", "interval1"))
    ale <- merge(ale, ale2, by = c("time", "interval2"))
    ale$ale <- ale$yhat_cumsum - ale$ale1 - ale$ale2 - ale$fJ0
    ale <- ale[order(ale$interval1, ale$interval2, ale$time), ]

    if (!center) {
        ale$ale <- ale$ale + mean_pred
    }

    interval_dists <- diff(quantile_vec1[c(1, 1:length(quantile_vec1), length(quantile_vec1))])
    interval_dists <- 0.5 * interval_dists

    ale$right <- quantile_vec1[ale$interval1 + 1] + interval_dists[ale$interval1 + 2]
    ale$left <- quantile_vec1[ale$interval1 + 1] - interval_dists[ale$interval1 + 1]

    interval_dists2 <- diff(quantile_vec2[c(1, 1:length(quantile_vec2), length(quantile_vec2))])
    interval_dists2 <- 0.5 * interval_dists2

    ale$bottom <- quantile_vec2[ale$interval2 + 1] + interval_dists2[ale$interval2 + 2]
    ale$top <- quantile_vec2[ale$interval2 + 1] - interval_dists2[ale$interval2 + 1]

    ale[, "_v1value_"] <- quantile_vec1[ale$interval1 + 1]
    ale[, "_v2value_"] <- quantile_vec2[ale$interval2 + 1]

    data.frame(
        "_v1name_" = var1,
        "_v2name_" = var2,
        "_v1type_" = "numerical",
        "_v2type_" = "numerical",
        "_v1value_" = ale$`_v1value_`,
        "_v2value_" = ale$`_v2value_`,
        "_times_" = ale$time,
        "_yhat_" = ale$ale,
        "_right_" = ale$right,
        "_left_" = ale$left,
        "_top_" = ale$top,
        "_bottom_" = ale$bottom,
        "_count_" = ifelse(is.na(ale$count), 0, ale$count),
        "_label_" = label,
        check.names = FALSE
    )
}
