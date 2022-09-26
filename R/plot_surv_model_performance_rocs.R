#' Plot ROC Curves for Survival Models
#'
#' This function plots objects of class `"surv_model_performance_rocs"` - ROC curves for
#' specific time points for survival models created using the `model_performance(..., type="roc")`.
#'
#' @param x an object of class `"surv_model_performance_rocs"` to be plotted
#' @param ... additional objects of class `"surv_model_performance_rocs"` to be plotted together
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, if `NULL` automaticaly generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' @param facet_ncol number of columns for arranging subplots
#'
#' @return An object of the class `ggplot`.
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
#' m_perf_roc <- model_performance(exp, type = "roc", times = c(100, 300))
#' plot(m_perf_roc)
#'
#' @export
plot.surv_model_performance_rocs <- function(x, ..., title = "ROC curves for selected timepoints", subtitle = NULL, colors = NULL, facet_ncol = NULL) {

    dfl <- c(list(x), list(...))

    alldfs <- lapply(dfl, function(x) {
        label <- attr(x, "label")

        cbind(x, label = label)
    })

    df <- do.call(rbind, alldfs)

    df$time_formatted <- paste0("t=", df$time)

    if (is.null(subtitle))
        subtitle <- paste0("created for the ", paste(unique(df$label), collapse = ", "), " model")

    num_colors <- length(unique(df$label))

    with(df, {
    ggplot(data = df, aes(x = FPR, y = TPR, group = label, color = label)) +
        geom_line(size = 0.8) +
        theme_drwhy() +
        xlab("1 - specificity (FPR)") +
        ylab("sensitivity (TPR)") +
        coord_fixed() +
        theme(panel.grid.major.x = element_line(color = "grey90", size = 0.5, linetype = 1),
        panel.grid.minor.x = element_line(color = "grey90", size = 0.5,  linetype = 1)) +
        labs(title = title, subtitle = subtitle) +
        scale_color_manual("", values = generate_discrete_color_scale(num_colors, colors)) +
        facet_wrap(~time, ncol = facet_ncol, labeller = function(x) lapply(x, function(x) paste0("t=", x)))
    })

}
