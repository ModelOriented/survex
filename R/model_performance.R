#' Dataset Level Performance Measures
#'
#' This function calculates metrics for survival models. The metrics calculated are C/D AUC, Brier score, and their integrated versions, as well as concordance index. It also can calculate ROC curves for specific selected time points.
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param ... other parameters, currently ignored
#' @param type character, either `"metrics"` or `"roc"`. If `"metrics"` then performance metrics are calculated, if `"roc"` ROC curves for selected time points are calculated.
#' @param times a numeric vector of times. If `type == "metrics"` then the survival function is evaluated at these times, if `type == "roc"` then the ROC curves are calculated at these times.
#'
#' @return An object of class `"model_performance_survival"`. It's a list of metric values calculated for the model. It contains:
#' - Harrell's concordance index [1]
#' - Brier score [2]
#' - C/D AUC using the estimator proposed by Uno et. al [3]
#' - integral of the Brier score
#' - integral of the C/D AUC
#'
#' @section References:
#' - [[1] Harrell Jr, Frank E., et al. "Regression modelling strategies for improved prognostic prediction." Statistics in medicine 3.2 (1984): 143-152.](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4780030207)
#' - [[2] Brier, Glenn W. "Verification of forecasts expressed in terms of probability." Monthly weather review 78.1 (1950): 1-3.](https://journals.ametsoc.org/view/journals/mwre/78/1/1520-0493_1950_078_0001_vofeit_2_0_co_2.xml)
#' - [[3] Uno, Hajime, et al. "Evaluating prediction rules for t-year survivors with censored regression models." Journal of the American Statistical Association 102.478 (2007): 527-537.](https://www.jstor.org/stable/27639883#metadata_info_tab_contents)
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' rsf_ranger <- ranger::ranger(Surv(time, status) ~ .,
#'                             data = veteran,
#'                             respect.unordered.factors = TRUE,
#'                             num.trees = 100,
#'                             mtry = 3,
#'                             max.depth = 5)
#'
#' rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ .,
#'                                 data = veteran)
#'
#' cph_exp <- explain(cph)
#' rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)],
#'                           y = Surv(veteran$time, veteran$status))
#' rsf_src_exp <- explain(rsf_src)
#'
#'
#' cph_model_performance <- model_performance(cph_exp)
#' rsf_ranger_model_performance <- model_performance(rsf_ranger_exp)
#' rsf_src_model_performance <- model_performance(rsf_src_exp)
#'
#' print(cph_model_performance)
#'
#' plot(rsf_ranger_model_performance, cph_model_performance,
#'      rsf_src_model_performance, metrics_type = "scalar")
#'
#' plot(rsf_ranger_model_performance, cph_model_performance, rsf_src_model_performance)
#'
#' cph_model_performance_roc <- model_performance(cph_exp, type = "roc", times = c(100, 500, 1200))
#' plot(cph_model_performance_roc)
#'
#' @rdname model_performance.surv_explainer
#' @export
model_performance <- function(explainer, ...) UseMethod("model_performance", explainer)

#' @rdname model_performance.surv_explainer
#' @export
model_performance.surv_explainer <- function(explainer,  ..., type = "metrics", times = NULL) {
    test_explainer(explainer, "model_performance", has_data = TRUE, has_y = TRUE, has_survival = TRUE, has_predict = TRUE)

    res <- surv_model_performance(explainer, ..., type = type, times = times)

    class(res) <- c("model_performance_survival", class(res))
    res
}

#' @export
model_performance.default <- DALEX::model_performance
