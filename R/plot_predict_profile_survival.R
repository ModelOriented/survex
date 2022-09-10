#' Plot Predict Profile for Survival Models
#'
#' This function plots objects of class `"predict_profile_survival"` - local explanations for survival models
#'
#' @param x an object of class  `"predict_profile_survival"` to be plotted
#' @param ... additional parameters passed to the `plot.surv_ceteris_paribus` function
#'
#' @return A grid of `ggplot` objects arranged with the `gridExtra::grid.arrange` function.
#'
#' @section Plot options:
#'
#' ## `plot.surv_ceteris_paribus`
#'
#' * `x` - an object of class `predict_profile_survival` to be plotted
#' * `...` - additional parameters, unused, currently ignored
#' * `colors` - character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' * `variable_type` - character, either `"numerical"`, `"categorical"` or `NULL` (default), select only one type of variable for plotting, or leave `NULL` for all
#' * `facet_ncol` - number of columns for arranging subplots
#' * `variables` - character, names of the variables to be plotted
#' * `numerical_plot_type` - character, either `"lines"`, or `"contours"` selects the type of numerical variable plots
#' * `title` - character, title of the plot
#' * `subtitle` - character, subtitle of the plot, if `NULL` automaticaly generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#'
#'
#' @family functions for plotting 'predict_profile_survival' objects
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' model <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' exp <- explain(model)
#'
#' p_profile <- predict_profile(exp, veteran[1, -c(3, 4)])
#'
#' plot(p_profile)
#'
#' p_profile_with_cat <- predict_profile(
#'     exp,
#'     veteran[1, -c(3, 4)],
#'     categorical_variables = c("trt", "prior")
#' )
#'
#' plot(p_profile_with_cat)
#' }
#' @export
plot.predict_profile_survival <- function(x, ...) {

    class(x) <- class(x)[-1]
    plot(x, ...)

}
