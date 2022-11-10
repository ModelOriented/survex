#' Plot Model Performance Metrics for Survival Models
#'
#' This function plots objects of class `"surv_model_performance"` - visualization of
#' metrics of different models created using the `model_performance(..., type="metrics")` function.
#'
#' @param x an object of class `"surv_model_performance"` to be plotted
#' @param ... additional objects of class `"surv_model_performance"` to be plotted together
#' @param metrics character, names of metrics to be plotted (subset of C/D AUC", "Brier score" for `metrics_type %in% c("time_dependent", "functional")` or subset of "C-index","Integrated Brier score", "Integrated C/D AUC" for `metrics_type == "scalar"`), by default (`NULL`) all metrics of a given type are plotted
#' @param metrics_type character, either one of `c("time_dependent","functional")` for functional metrics or `"scalar"` for scalar metrics
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param facet_ncol number of columns for arranging subplots
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
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
#' m_perf <- model_performance(exp)
#' plot(m_perf)
#'
#' @export
plot.surv_model_performance <- function(x,
                                        ...,
                                        metrics = NULL,
                                        metrics_type = "time_dependent",
                                        title = "Model performance",
                                        subtitle = "default",
                                        facet_ncol = NULL,
                                        colors = NULL) {


    if (metrics_type %in% c("time_dependent", "functional")) {
        pl <- plot_td_surv_model_performance(x, ..., metrics = metrics, title = title, subtitle = subtitle, facet_ncol = facet_ncol, colors = colors)
    }


    if (metrics_type == "scalar") {
        pl <- plot_scalar_surv_model_performance(x, ..., metrics = metrics, title = title, subtitle = subtitle, facet_ncol = facet_ncol, colors = colors)
    }

    pl

}


#' @importFrom DALEX theme_drwhy
plot_td_surv_model_performance <- function(x, ..., metrics = NULL, title = NULL, subtitle = "default", facet_ncol = NULL, colors = NULL) {

    df <- concatenate_td_dfs(x, ...)

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", paste(unique(df$label), collapse = ", "), " model")
    }

    if (is.null(metrics)) metrics <- c("C/D AUC", "Brier score")

    num_colors <- length(unique(df$label))

    ggplot(data = df[df$ind %in% metrics, ], aes_string(x = "times", y = "values", group = "label", color = "label")) +
        geom_line(size = 0.8) +
        theme_drwhy() +
        xlab("") +
        ylab("metric value") +
        labs(title = title, subtitle = subtitle) +
        scale_color_manual("", values = generate_discrete_color_scale(num_colors, colors)) +
        facet_wrap(~ind, ncol = facet_ncol, scales = "free_y")

}

#' @importFrom DALEX theme_drwhy
plot_scalar_surv_model_performance <- function(x, ..., metrics = NULL, title = NULL, subtitle = NULL, facet_ncol = NULL, colors = NULL) {
    df <- concatenate_dfs(x, ...)

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", paste(unique(df$label), collapse = ", "), " model")
    }

    if (!is.null(metrics)) df <- df[df$ind %in% metrics, ]

    num_colors <- length(unique(df$label))

    ggplot(data = df, aes_string(x = "label", y = "values", fill = "label")) +
        geom_col() +
        theme_drwhy() +
        xlab("") +
        ylab("metric value") +
        labs(title = title, subtitle = subtitle) +
        scale_fill_manual("", values = generate_discrete_color_scale(num_colors, colors)) +
        facet_wrap(~ind, ncol = facet_ncol, scales = "free_y")


}


concatenate_td_dfs <- function(x, ...) {
    all_things <- c(list(x), list(...))

    all_dfs <- lapply(all_things, function(x) {
        df <- data.frame(`Brier score` = x$brier_score,
                         `C/D AUC` = x$auc, check.names = FALSE)

        df <- stack(df)
        times <-  rep(x$eval_times, 2)
        label <-  attr(x, "label")
        df <- cbind(times, df, label)
    })

    do.call(rbind, all_dfs)

}


concatenate_dfs <- function(x, ...) {
    all_things <- c(list(x), list(...))

    all_dfs <- lapply(all_things, function(x) {
        tmp_list <- lapply(x, function(metric) {
            if(!is.null(attr(metric, "loss_type"))){
               if(attr(metric, "loss_type") != "time-dependent"){
                metric[1]}
            }
            })
        tmp_list[sapply(tmp_list, is.null)] <- NULL
        df <- data.frame(tmp_list,
                         check.names = FALSE)
        df <- stack(df)
        label <- attr(x, "label")
        df <- cbind(df, label)
    })

    do.call(rbind, all_dfs)
}
