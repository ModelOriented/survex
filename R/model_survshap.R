#'@export
model_survshap <-
    function(explainer, ...)
        UseMethod("model_survshap", explainer)

#'@export
model_survshap.surv_explainer <- function(explainer,
                                          calculation_method = "kernelshap",
                                          aggregation_method = "integral",
                                          ...,
                                          N = NULL) {

    test_explainer(
        explainer,
        has_data = TRUE,
        has_y = TRUE,
        has_survival = TRUE,
        function_name = "model_survshap"
    )
    observations <- explainer$data

    if (!is.null(N)) {
        selected_observations <- sample(1:nrow(observations), N)
        observations <- selected_observations
    }

    shap_values <- surv_shap(
        explainer = explainer,
        new_observation = observations,
        calculation_method = calculation_method,
        aggregation_method = aggregation_method
    )

    attr(shap_values, "label") <- explainer$label
    shap_values$event_times <- explainer$y[explainer$y[, 1] <= max(explainer$times), 1]
    shap_values$event_statuses <- explainer$y[explainer$y[, 1] <= max(explainer$times), 2]
    return(shap_values)

}
