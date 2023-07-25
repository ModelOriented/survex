#' Plot SurvSHAP(t) Explanations for Survival Models
#'
#' This functions plots objects of class `surv_shap` - time-dependent explanations of
#' survival models created using the `predict_parts(..., type="survshap")` function.
#'
#' @param x an object of class `"surv_shap"` to be plotted
#' @param ... additional objects of class `surv_shap` to be plotted together
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param max_vars maximum number of variables to be plotted (least important variables are ignored)
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' @param rug character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`.
#' @param rug_colors character vector containing two colors (containing either hex codes "#FF69B4", or names "blue"). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
#'
#' @return An object of the class `ggplot`.
#'
#' @family functions for plotting 'predict_parts_survival' objects
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' model <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' exp <- explain(model)
#'
#' p_parts_shap <- predict_parts(exp, veteran[1, -c(3, 4)], type = "survshap")
#' plot(p_parts_shap)
#' }
#'
#' @export
plot.surv_shap <- function(x,
                           ...,
                           title = "SurvSHAP(t)",
                           subtitle = "default",
                           max_vars = 7,
                           colors = NULL,
                           rug = "all",
                           rug_colors = c("#dd0000", "#222222")) {

    dfl <- c(list(x), list(...))

    long_df <- lapply(dfl, function(x) {
        label <- attr(x, "label")
        cols <- sort(head(order(x$aggregate, decreasing = TRUE), max_vars))
        sv <- x$result[,cols]
        times <- x$eval_times
        transposed <- as.data.frame(cbind(times = times, sv))
        rownames(transposed) <- NULL
        long_df <- cbind(
            times = transposed$times,
            stack(transposed, select = -times),
            label = label
        )
    })

    transformed_rug_dfs <- lapply(dfl, function(x){
        label <- attr(x, "label")
        rug_df <- data.frame(times = x$event_times, statuses = as.character(x$event_statuses), label = label)
    })

    rug_df <- do.call(rbind, transformed_rug_dfs)

    long_df <- do.call(rbind, long_df)
    label <- unique(long_df$label)

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", paste(label, collapse = ", "), " model")
    }

    n_colors <- length(unique(long_df$ind))

    y_lab <- "SurvSHAP(t) value"


    base_plot <- with(long_df, {
    ggplot(data = long_df, aes(x = times, y = values, color = ind)) +
        geom_line(linewidth = 0.8, size = 0.8) +
        ylab(y_lab) + xlab("") +
        xlim(c(0,NA))+
        labs(title = title, subtitle = subtitle) +
        scale_color_manual("variable", values = generate_discrete_color_scale(n_colors, colors)) +
        theme_default_survex() +
        facet_wrap(~label, ncol = 1, scales = "free_y")
    })

    return_plot <- add_rug_to_plot(base_plot, rug_df, rug, rug_colors)

    return(return_plot)

}



#'@export
plot.aggregated_surv_shap <- function(x,
                                      ...,
                                      title = "Feature importance according to aggregated |SurvSHAP(t)|",
                                      subtitle = "default",
                                      xlab_left = "Importance",
                                      ylab_right = "Aggregated SurvSHAP(t) value",
                                      max_vars = 7,
                                      colors = NULL,
                                      rug = "all",
                                      rug_colors = c("#dd0000", "#222222")){


    old_x <- x
    x$result <- aggregate_shap_multiple_observations(x$result, colnames(x$result[[1]]), function(x) mean(abs(x)))
    x$aggregate <- apply(do.call(rbind, x$aggregate), 2, function(x) mean(abs(x)))

    right_plot <- plot.surv_shap(x = x,
                                 ... = ...,
                                 title = NULL,
                                 subtitle = NULL,
                                 max_vars = max_vars,
                                 colors = colors,
                                 rug = rug,
                                 rug_colors = rug_colors) +
                  labs(y = ylab_right)


    dfl <- c(list(x), list(...))

    df_list <- lapply(dfl, function(x) {
        label <- attr(x, "label")
        values <- x$aggregate
        vars <- names(x$aggregate)
        df <- data.frame(label, values, vars)
        rownames(df) <- NULL
        df
    })

    long_df <- do.call("rbind", df_list)
    long_df <- long_df[order(long_df$values, decreasing = TRUE),]

    label <- unique(long_df$label)
    subtitle <- paste0(
        "created for the ", paste(label, collapse = ", "), " model ",
        "(n=", x$n_observations, ")"
    )

    left_plot <- with(long_df, {
        ggplot(long_df, aes(x = values, y = reorder(vars, values))) +
            geom_col(fill = "#46bac2") +
            theme_default_survex() +
            facet_wrap(~label, ncol = 1, scales = "free_y") +
            labs(x = xlab_left) +
            theme(axis.title.y = element_blank())

    })


    pl <- left_plot +
          right_plot +
          patchwork::plot_layout(widths = c(3,5), guides = "collect") +
          patchwork::plot_annotation(title = title, subtitle = subtitle) &
          theme(legend.position = "top",
                plot.title = element_text(color = "#371ea3", size = 16, hjust = 0),
                plot.subtitle = element_text(color = "#371ea3", hjust = 0),)

    pl
}
