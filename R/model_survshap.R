#' Global SHAP Values
#'
#' This function computes global SHAP values.
#'
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
#' @examples
#' \donttest{
#' veteran <- survival::veteran
#' rsf_ranger <- ranger::ranger(
#'     survival::Surv(time, status) ~ .,
#'     data = veteran,
#'     respect.unordered.factors = TRUE,
#'     num.trees = 100,
#'     mtry = 3,
#'     max.depth = 5
#' )
#' rsf_ranger_exp <- explain(
#'     rsf_ranger,
#'     data = veteran[, -c(3, 4)],
#'     y = survival::Surv(veteran$time, veteran$status),
#'     verbose = FALSE
#' )
#'
#' ranger_global_survshap <- model_survshap(
#'     explainer = rsf_ranger_exp,
#'     new_observation = veteran[
#'         c(1:4, 17:20, 110:113, 126:129),
#'         !colnames(veteran) %in% c("time", "status")
#'     ],
#'     y_true = survival::Surv(
#'         veteran$time[c(1:4, 17:20, 110:113, 126:129)],
#'         veteran$status[c(1:4, 17:20, 110:113, 126:129)]
#'     ),
#'     aggregation_method = "integral",
#'     calculation_method = "kernelshap",
#' )
#' plot(ranger_global_survshap)
#' plot(ranger_global_survshap, geom = "beeswarm")
#' plot(ranger_global_survshap, geom = "profile", color_variable = "karno")
#' }
#'
#' @rdname model_survshap.surv_explainer
#' @export
model_survshap <- function(explainer, ...) {
        UseMethod("model_survshap", explainer)
    }

#' @rdname model_survshap.surv_explainer
#' @export
model_survshap.surv_explainer <- function(explainer,
                                          new_observation = NULL,
                                          y_true = NULL,
                                          N = NULL,
                                          calculation_method = "kernelshap",
                                          aggregation_method = "integral",
                                          output_type = "survival",
                                          ...) {
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
        y_true <- explainer$y
    }

    shap_values <- surv_shap(
        explainer = explainer,
        new_observation = observations,
        output_type = output_type,
        N = N,
        y_true = y_true,
        calculation_method = calculation_method,
        aggregation_method = aggregation_method
    )

    attr(shap_values, "label") <- explainer$label
    shap_values$event_times <- explainer$y[explainer$y[, 1] <= max(explainer$times), 1]
    shap_values$event_statuses <- explainer$y[explainer$y[, 1] <= max(explainer$times), 2]
    shap_values$output_type <- output_type
    return(shap_values)
}
