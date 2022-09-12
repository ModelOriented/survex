#' Model Predictions for Survival Models
#'
#' This function allows for calculating model prediction in a unified way.
#'
#' @param object a model to make the predictions, preprocessed by the `explain()` function.
#' @param newdata data used for the prediction
#' @param output_type character, either `"risk"`, `"survival"` or `"chf"` depending on the desired output
#' @param times a numeric vector of times for the survival and cumulative hazard function predictions to be evaluated at. If `"output_type == "risk"` this argument is ignored, if left `NULL` then it is extracted from `object$times`.
#' @param ... other arguments, currently ignored
#'
#' @return A vector or matrix containing the prediction.
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' rsf_ranger <- ranger::ranger(Surv(time, status) ~ .,
#'                              data = veteran,
#'                              respect.unordered.factors = TRUE,
#'                              num.trees = 100,
#'                              mtry = 3,
#'                              max.depth = 5)
#'
#' cph_exp <- explain(cph)
#'
#' rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)],
#'                           y = Surv(veteran$time, veteran$status))
#'
#'
#' predict(cph_exp, veteran[1, ], output_type = "survival")[, 1:10]
#'
#' predict(cph_exp, veteran[1, ], output_type = "risk")
#'
#' predict(rsf_ranger_exp, veteran[1, ], output_type = "chf")[, 1:10]
#'
#' @export
predict.surv_explainer <- function(object, newdata = NULL, output_type = "survival", times = NULL,  ...) {

    if (is.null(newdata)) newdata <- object$data
    if (is.null(times)) times <- object$times

    model <- object$model


    switch(output_type,
        "risk" = object$predict_function(model, newdata),
        "survival" = object$predict_survival_function(model, newdata, times),
        "chf" = object$predict_cumulative_hazard_function(model, newdata, times),
        stop("`output_type` should be one of `risk`, `survival` or `chf`")
    )
}
