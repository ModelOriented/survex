#' Plot Permuatational Feature Importance for Survival Models
#'
#' This function plots feature importance objects created for survival models.
#'
#' @param x an object of class `"surv_feature_importance"` to be plotted
#' @param ... additional objects of class `"surv_feature_importance"` to be plotted together
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, if `NULL` automaticaly generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#' @param max_vars maximum number of variables to be plotted (least important variables are ignored)
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#'
#' @return A `ggplot2` plot.
#'
#' @family functions for plotting 'model_parts_survival' objects
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' model <- coxph(Surv(time, status) ~ ., data = veteran, x = TRUE, model = TRUE, y = TRUE)
#' model_rf <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' explainer <- explain(model)
#' explainer_rf <- explain(model_rf)
#'
#' mp <- model_parts(explainer)
#' mp_rf <- model_parts(explainer_rf)
#'
#' plot(mp, mp_rf)
#' }
#' @importFrom utils stack head
#' @importFrom DALEX theme_drwhy
#' @export
plot.surv_feature_importance <- function(x, ...,
                                                    title = "Time-dependent feature importance",
                                                    subtitle = NULL,
                                                    max_vars = 6,
                                                    colors = NULL) {

    df_list <- c(list(x), list(...))

    transformed_dfs <- lapply(df_list, function(x) {
        x <- x$result
        label <- unique(x$label)
        x <- x[x$permutation == 0, !colnames(x) %in% c("permutation", "label", "_baseline_")]
        plotting_df <- with(x, cbind(x[1], stack(x, select = -times), label, row.names = NULL))
    })


    plotting_df <- do.call(rbind, transformed_dfs)

    label <- unique(plotting_df$label)

    subs <- aggregate(plotting_df$value, by = list(var = plotting_df$ind), function(x) sum(abs(x)))

    subs <- subs[order(subs$x, decreasing = TRUE), ]
    plotting_df <- plotting_df[plotting_df$ind %in% c("_full_model_", as.character(head(subs$var, max_vars))), ]

    num_variables <- length(unique(plotting_df$ind))

    additional_info <- switch(attr(x, "type"),
                              "raw" = "",
                              "ratio" = "\ndivided by the loss of full model",
                              "difference" = "\nwith loss of full model subtracted")

    if (!is.null(attr(x, "loss_name"))) {
        y_lab <- paste0(paste(attr(x, "loss_name")[1], "loss after permutations"), additional_info)
    } else {
        y_lab <- paste0("Loss function after variable's permutations", additional_info)
    }

    if (is.null(subtitle)) {
        glm_labels <- paste0(label, collapse = ", ")
        subtitle <- paste0("created for the ", glm_labels, " model")
    }


    ggplot(data = plotting_df, aes_string(x = "times", y = "values", color = "ind", label = "ind")) +
        geom_line(size = 0.8) +
        theme_drwhy() +
        xlab("") +
        ylab(y_lab) +
        scale_color_manual(name = "Variable", values = c("#000000", generate_discrete_color_scale(num_variables, colors))) +
        labs(title = title, subtitle = subtitle) +
        facet_wrap(~label)

}
