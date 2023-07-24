#' Global SHAP Values
#'
#' This function computes global SHAP values.
#'
#' @param N A positive integer indicating the number of observations that should be used to compute global SHAP values.
#' @inheritParams surv_shap
#'
#' @details
#' If specifying `y_true`, also `new_observation` must be specified.
#' Using the argument `new_observation`, global SHAP values are computed for the provided data. Otherwise,
#' global SHAP values are computed for the data, the `explainer` was trained with.
#'
#'
#' @return An object of class `aggregated_surv_shap` containing the computed global SHAP values.
#'
#' @rdname model_survshap.surv_explainer
#' @export
model_survshap <-
    function(explainer, ...)
        UseMethod("model_survshap", explainer)

#' @rdname model_survshap.surv_explainer
#' @export
model_survshap.surv_explainer <- function(explainer,
                                          calculation_method = "kernelshap",
                                          aggregation_method = "integral",
                                          new_observation = NULL,
                                          y_true = NULL,
                                          ...,
                                          N = NULL) {

    stopifnot(
        "`N` must be a positive integer" = ifelse(
            !is.null(N),
            is.integer(N) && N > 0L,
            TRUE
        ),
        "`y_true` must be either a matrix with one per observation in `new_observation` or a vector of length == 2" = ifelse(
            !is.null(y_true),
            ifelse(
                is.matrix(y_true),
                nrow(new_observation) == nrow(y_true),
                is.null(dim(y_true)) && length(y_true) == 2L
            ),
            TRUE
        ))

    test_explainer(
        explainer,
        has_data = TRUE,
        has_y = TRUE,
        has_survival = TRUE,
        function_name = "model_survshap"
    )
    if (!is.null(new_observation)) {
        if (length(setdiff(colnames(new_observation), colnames(explainer$data))) == 0) {
            observations <- new_observation
            y_true <- y_true
        } else {
            stop("`new_observation` must have the same column names as the data, the `explainer` was trained with.")
        }
    } else {
        observations <- explainer$data
        y_true <- NULL
    }

    if (!is.null(N)) {
        selected_observations <- sample(1:nrow(observations), N)
        observations <- observations[selected_observations, ]
    }

    shap_values <- surv_shap(
        explainer = explainer,
        new_observation = observations,
        y_true = y_true,
        calculation_method = calculation_method,
        aggregation_method = aggregation_method
    )

    attr(shap_values, "label") <- explainer$label
    shap_values$event_times <- explainer$y[explainer$y[, 1] <= max(explainer$times), 1]
    shap_values$event_statuses <- explainer$y[explainer$y[, 1] <= max(explainer$times), 2]
    return(shap_values)

}
