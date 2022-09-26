#' Plot Model Performance for Survival Models
#'
#' This function is a wrapper for plotting `model_performance` objects created for survival models
#' using the `model_performance()` function.
#'
#' @param x an object of class `"model_performance_survival"` to be plotted
#' @param ... additional parameters passed to the `plot.surv_model_performance` or `plot.surv_model_performance_rocs` function
#'
#' @return An object of the class `ggplot`.
#'
#' @section Plot options:
#'
#' ## `plot.surv_model_performance`
#'
#' * `x` - an object of class `"surv_model_performance"` to be plotted
#' * `...` - additional objects of class `"surv_model_performance"` to be plotted together
#' * `metrics` - character, names of metrics to be plotted (subset of C/D AUC", "Brier score" for `metrics_type %in% c("time_dependent", "functional")` or subset of "C-index","Integrated Brier score", "Integrated C/D AUC" for `metrics_type == "scalar"`), by default (`NULL`) all metrics of a given type are plotted
#' * `metrics_type` - character, either one of `c("time_dependent","functional")` for functional metrics or `"scalar"` for scalar metrics
#' * `title` - character, title of the plot
#' * `subtitle` - character, subtitle of the plot, if `NULL` automaticaly generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#' * `facet_ncol` - number of columns for arranging subplots
#' * `colors` - character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#'
#' ## `plot.surv_model_performance_rocs`
#'
#' * `x` - an object of class `"surv_model_performance_rocs"` to be plotted
#' * `...` - additional objects of class `"surv_model_performance_rocs"` to be plotted together
#' * `title` - character, title of the plot
#' * `subtitle` - character, subtitle of the plot, if `NULL` automaticaly generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#' * `colors` - character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' * `facet_ncol` - number of columns for arranging subplots
#'
#' @family functions for plotting 'model_performance_survival' objects
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#' model <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' exp <- explain(model)
#'
#' m_perf <- model_performance(exp)
#' plot(m_perf, metrics_type = "functional")
#'
#' m_perf_roc <- model_performance(exp, type = "roc", times = c(100, 300))
#' plot(m_perf_roc)
#'
#' @export
plot.model_performance_survival <- function(x, ...) {
    class(x) <- class(x)[-1]
    plot(x, ...)
}
