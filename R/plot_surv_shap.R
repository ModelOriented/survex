#' Plot SurvSHAP(t) Explanations for Survival Models
#'
#' This functions plots objects of class `surv_shap` - time-dependent explanations of
#' survival models created using the `predict_parts(..., type="survshap")` function.
#'
#' @param x an object of class `"surv_shap"` to be plotted
#' @param ... additional objects of class `surv_shap` to be plotted together
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' @param rug character vector. One of "all", "events", "censors", "none" or NULL. Which times to mark on the x axis in `geom_rug()`.
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
                           colors = NULL,
                           rug = "all",
                           rug_colors = c("#dd0000", "#222222")) {

    dfl <- c(list(x), list(...))

    long_df <- lapply(dfl, function(x) {
        label <- attr(x, "label")
        sv <- x$result
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
        labs(title = title, subtitle = subtitle) +
        scale_color_manual("variable", values = generate_discrete_color_scale(n_colors, colors)) +
        theme_drwhy() +
        facet_wrap(~label, ncol = 1, scales = "free_y")
    })

    if (rug == "all"){
        return_plot <- base_plot +
            geom_rug(data = rug_df[rug_df$statuses == 1,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[1]) +
            geom_rug(data = rug_df[rug_df$statuses == 0,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[2])
    } else if (rug == "events") {
        return_plot <- base_plot +
            geom_rug(data = rug_df[rug_df$statuses == 1,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[1])
    } else if (rug == "censors") {
        return_plot <- base_plot +
            geom_rug(data = rug_df[rug_df$statuses == 0,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[2])
    } else {
        return_plot <- base_plot
    }
    return(return_plot)

}
