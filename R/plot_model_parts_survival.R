#' Plot Model Parts for Survival Models
#'
#' This function is a wrapper for plotting `model_parts` objects created for survival models
#' using the `model_parts()` function.
#'
#' @param x an object of class `"model_parts_survival"` to be plotted
#' @param ... additional parameters passed to the `plot.surv_feature_importance` function
#'
#' @return An object of the class `ggplot`.
#'
#' @section Plot options:
#'
#' * `title` - character, title of the plot
#' * `subtitle` - character, subtitle of the plot, if `NULL` automatically generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#' * `max_vars` - maximum number of variables to be plotted (least important variables are ignored)
#' * `colors` - character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' * `rug` - character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`.
#' * `rug_colors` - character vector containing two colors (containing either hex codes "#FF69B4", or names "blue"). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
#'
#' @family functions for plotting 'model_parts_survival' objects
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#' model <- coxph(Surv(time, status) ~ ., data = veteran, x = TRUE, model = TRUE, y = TRUE)
#' explainer <- explain(model)
#'
#' mp <- model_parts(explainer)
#'
#' plot(mp)
#'
#' @export
plot.model_parts_survival <- function(x, ...) {
    class(x) <- class(x)[-1]
    plot(x, ...)
}
