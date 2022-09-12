#' Helper functions for `predict_parts.R`
#'
#' @param explainer a model to be explained, preprocessed by the `explain` function
#' @param new_observation a new observation for which predictions need to be explained
#' @param ... additional parameters, passed to internal functions
#' @param y_true a two element numeric vector or matrix of one row and two columns, the first element being the true observed time and the second the status of the observation, used for plotting
#' @param calculation_method a character, only `"kernel"` is implemented for now.
#' @param aggregation_method a character, either `"mean_absolute"` or `"integral"`, `"max_absolute"`, `"sum_of_squares"`
#' @param path ignored, placeholder the not implemented `"sampling"` method
#' @param B ignored, placeholder the not implemented `"sampling"` method
#' @param exact ignored, placeholder the not implemented `"sampling"` method
#'
#' @return A list, containing the calculated SurvSHAP(t) results in the `result` field
#'
#'
#' @keywords internal
surv_shap <- function(explainer,
                      new_observation,
                      ...,
                      y_true = NULL,

                      calculation_method = "kernel",
                      aggregation_method = "integral",
                      path = "average",
                      B = 25,
                      exact = FALSE
) {
    test_explainer(explainer, "surv_shap", has_data = TRUE, has_y = TRUE, has_survival = TRUE)
    new_observation <- new_observation[, !colnames(new_observation) %in% colnames(explainer$y)]
    if (ncol(explainer$data) != ncol(new_observation)) stop("New observation and data have different number of columns(variables)")

    event_inds <- explainer$y[, 2]
    event_times <- explainer$y[, 1]

    if (!is.null(y_true)) {
        if (is.matrix(y_true)) {
            y_true_ind <- y_true[1, 2]
            y_true_time <- y_true[1, 1]
        } else {
            y_true_ind <- y_true[2]
            y_true_time <- y_true[1]
        }
    }

    res <- switch(calculation_method,
                  "kernel" = shap_kernel(explainer, new_observation, aggregation_method, ...),
                  stop("Only calculation method = `kernel` is implemented"))

    if (!is.null(y_true)) res$y_true <- c(y_true_time = y_true_time, y_true_ind = y_true_ind)

    res$aggregate <- aggregate_surv_shap(res, aggregation_method)

    class(res) <- "surv_shap"
    res
}


shap_kernel <- function(explainer, new_observation, aggregation_method,  ...) {


    timestamps <- explainer$times


    p <- ncol(explainer$data)

    target_sf <- explainer$predict_survival_function(explainer$model, new_observation, timestamps)
    sfs <- explainer$predict_survival_function(explainer$model, explainer$data, timestamps)
    baseline_sf <- apply(sfs, 2, mean)


    permutations <- expand.grid(rep(list(0:1), p))
    kernel_weights <- generate_shap_kernel_weights(permutations, p)

    shap_values <- calculate_shap_values(explainer, explainer$model, baseline_sf, explainer$data, permutations, kernel_weights, new_observation, timestamps)

    shap_values <- as.data.frame(shap_values, row.names = colnames(explainer$data))
    colnames(shap_values) <- paste("t=", timestamps, sep = "")


    ret <- list()
    ret$eval_times <- timestamps
    ret$result <- t(shap_values)
    ret$variable_values <- new_observation

    return(ret)
}

generate_shap_kernel_weights <- function(permutations, p) {

    apply(permutations, 1, function(row) {
        row <- as.numeric(row)
        num_available_variables = sum(row != 0)

        if (num_available_variables == 0 || num_available_variables == p) 1e12
        else {
            (p - 1) / (choose(p, num_available_variables) * num_available_variables * (p - num_available_variables))
        }
        })
}


calculate_shap_values <- function(explainer, model, avg_survival_function, data, simplified_inputs, shap_kernel_weights, new_observation, timestamps) {

    W <- diag(shap_kernel_weights)
    X <- as.matrix(simplified_inputs)


    R <- solve((t(X) %*% W %*% X)) %*% (t(X) %*% W)

    y <- make_prediction_for_simplified_input(explainer, model, data, simplified_inputs, new_observation, timestamps)

    y <- sweep(y,
               2,
               avg_survival_function)

    R %*% y

}

make_prediction_for_simplified_input <- function(explainer, model, data, simplified_inputs, new_observation, timestamps) {


    preds <- apply(simplified_inputs, 1, function(row) {
        row <- as.logical(row)

        X_tmp <- data
        X_tmp[, row] <- new_observation[, row]
        colnames(X_tmp) <- colnames(data)

        apply(explainer$predict_survival_function(model, X_tmp, timestamps), 2, mean)

    })

    return(t(preds))


}

aggregate_surv_shap <- function(survshap, method) {

    switch(method,
           "sum_of_squares" = return(apply(survshap$result, 2, function(x) sum(x^2))),
           "mean_absolute" = return(apply(survshap$result, 2, function(x) mean(abs(x)))),
           "max_absolute" = return(apply(survshap$result, 2, function(x) max(abs(x)))),
           "integral" = return(apply(survshap$result, 2, function(x) {
               x <- abs(x)
               names(x) <- NULL
               times <- survshap$eval_times
               n <- length(x)
               i <- (x[1:(n - 1)] + x[2:n]) * diff(times) / 2
               cumsum(c(0, i))[length(cumsum(c(0, i)))] / (max(times) - min(times))
           })),
           stop("aggregation_method has to be one of `sum_of_squares`, `mean_absolute`, `max_absolute` or `integral`"))

}
