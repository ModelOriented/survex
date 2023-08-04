#' Plot SurvSHAP(t) Explanations for Survival Models
#'
#' This functions plots objects of class `surv_shap` - time-dependent explanations of
#' survival models created using the `predict_parts(..., type="survshap")` function.
#'
#' @param x an object of class `surv_shap` to be plotted
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
        geom_line(linewidth = 0.8) +
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


#' Plot Aggregated SurvSHAP(t) Explanations for Survival Models
#'
#' This functions plots objects of class `aggregated_surv_shap` - aggregated time-dependent
#' explanations of survival models created using the `model_survshap()` function.
#'
#' @param x an object of class `aggregated_surv_shap` to be plotted
#' @param kind character, one of `"importance"`, `"swarm"`, or `"profile"`. Type of chart to be plotted: `"importance"` shows the importance of variables over time and aggregated, `"swarm"` shows the distribution of SurvSHAP(t) values for variables and observations, `"profile"` shows the dependence of SurvSHAP(t) values on variable values.
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for the XXX model (n = YYY)", where XXX is the explainer label and YYY is the number of observations used for calculations
#' @param max_vars maximum number of variables to be plotted (least important variables are ignored), by default 7
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#'
#' @return An object of the class `ggplot`.
#'
#' @section Plot options:
#'
#' ## `plot.aggregated_surv_shap(type = "importance")`
#'
#' * `rug` - character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`.
#' * `rug_colors` - character vector containing two colors (containing either hex codes "#FF69B4", or names "blue"). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
#' * `xlab_left, ylab_right` - axis labels for left and right plots (due to different aggregation possibilities)
#'
#'
#' ## `plot.aggregated_surv_shap(type = "swarm")`
#'
#' * no additional options
#'
#'
#' ## `plot.aggregated_surv_shap(type = "swarm")`
#'
#' * `variable` - variable for which the profile is to be plotted, by default first from result data
#' * `color_variable` - variable used to denote the color, by default equal to `variable`
#'
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
#'@export
plot.aggregated_surv_shap <- function(x,
                                      kind = "importance",
                                      ...,
                                      colors = NULL){
    if (is.null(colors)){
        colors <- c(low = "#9fe5bd",
                    mid = "#46bac2",
                    high = "#371ea3")
    }

    switch(
        kind,
        "importance" = plot_shap_global_importance(x = x,
                                          ... = ...,
                                          colors = colors),
        "swarm" = plot_shap_global_swarm(x = x,
                                     ... = ...,
                                     colors = colors),
        "profile" = plot_shap_global_profile(x = x,
                                           ... = ...,
                                           colors = colors),
        stop("`kind` must be one of 'importance', 'swarm' or 'profile'")
    )
}

plot_shap_global_importance <- function(x,
                                        ...,
                                        title = "Feature importance according to aggregated |SurvSHAP(t)|",
                                        subtitle = "default",
                                        max_vars = 7,
                                        colors = NULL,
                                        rug = "all",
                                        rug_colors = c("#dd0000", "#222222"),
                                        xlab_left = "Average |aggregated SurvSHAP(t)| value",
                                        ylab_right = "Average |SurvSHAP(t)| value"){

    x$result <- aggregate_shap_multiple_observations(x$result, colnames(x$result[[1]]), function(x) mean(abs(x)))
    x$aggregate <- apply(do.call(rbind, x$aggregate), 2, function(x) mean(abs(x)))

    right_plot <- plot.surv_shap(x = x,
                                 title = NULL,
                                 subtitle = NULL,
                                 max_vars = max_vars,
                                 colors = NULL,
                                 rug = rug,
                                 rug_colors = rug_colors) +
        labs(y = ylab_right)

    label <- attr(x, "label")
    long_df <- stack(x$aggregate)
    long_df <- long_df[order(long_df$values, decreasing = TRUE),][1:min(max_vars, length(x$aggregate)), ]

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0(
            "created for the ", label, " model ",
            "(n=", x$n_observations, ")"
        )
    }

    left_plot <- with(long_df, {
        ggplot(long_df, aes(x = values, y = reorder(ind, values))) +
            geom_col(fill = colors[2]) +
            theme_default_survex() +
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

    return(pl)
}

plot_shap_global_swarm <- function(x,
                                   ...,
                                   title = "Aggregated SurvSHAP(t) values summary",
                                   subtitle = "default",
                                   max_vars = 7,
                                   colors = NULL){

    df <- as.data.frame(do.call(rbind, x$aggregate))
    cols <- names(sort(colMeans(abs(df))))[1:min(max_vars, length(df))]
    df <- df[,cols]
    df <- stack(df)
    colnames(df) <- c("shap_value", "variable")

    original_values <- as.data.frame(x$variable_values)[,cols]
    var_value <- preprocess_values_to_common_scale(original_values)
    df <- cbind(df, var_value)

    label <- attr(x, "label")
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0(
            "created for the ", label, " model ",
            "(n=", x$n_observations, ")"
        )
    }
    with(df, {
    ggplot(data = df, aes(x = shap_value, y = variable, color = var_value)) +
        geom_vline(xintercept = 0, color = "#ceced9", linetype="solid") +
        geom_jitter(width=0) +
        scale_color_gradient2(
            name = "Variable value",
            low = colors[1],
            mid = colors[2],
            high = colors[3],
            midpoint = 0.5,
            limits=c(0,1),
            breaks = c(0, 1),
            labels=c("", "")) +
        labs(title = title, subtitle = subtitle,
             x = "Aggregated SurvSHAP(t) value",
             y = "Variable") +
        theme_default_survex() +
        theme(legend.position = "bottom") +
        guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))
    }
  )
}

plot_shap_global_profile <- function(x,
                                     ...,
                                     variable = NULL,
                                     color_variable = NULL,
                                     title = "Aggregated SurvSHAP(t) profile",
                                     subtitle = "default",
                                     max_vars = 7,
                                     colors = NULL){

    df <- as.data.frame(do.call(rbind, x$aggregate))

    if (is.null(variable)){
        variable <- colnames(df)[1]
    }
    if (is.null(color_variable)){
        color_variable <- variable
    }

    shap_val <- df[,variable]

    original_values <- as.data.frame(x$variable_values)
    var_vals <- original_values[,c(variable, color_variable)]

    df <- cbind(shap_val, var_vals)
    colnames(df) <- c("shap_val", "variable_val", "color_variable_val")

    label <- attr(x, "label")
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0(
            "created for the ", label, " model ",
            "(n=", x$n_observations, ")"
        )
    }

    p <- with(df, {
      ggplot(df, aes(x = variable_val, y = shap_val, color = color_variable_val)) +
        geom_hline(yintercept = 0, color = "#ceced9", linetype="solid") +
        geom_point() +
        geom_rug(aes(x = df$variable_val), inherit.aes=F, color = "#ceced9") +
        labs(x = paste(variable, "value"),
             y = "Aggregated SurvSHAP(t) value",
             title = title,
             subtitle = subtitle) +
        theme_default_survex() +
        theme(legend.position = "bottom")
    })

    if (!is.factor(df$color_variable_val)) {
        p + scale_color_gradient2(
            name = paste(color_variable, "value"),
            low = colors[1],
            mid = colors[2],
            high = colors[3],
            midpoint = median(df$color_variable_val))
    } else {
        p + scale_color_manual(name = paste(color_variable, "value"),
                           values = generate_discrete_color_scale(length(unique(df$color_variable_val)), colors))
    }
}

preprocess_values_to_common_scale <- function(data) {
    # Scale numerical columns to range [0, 1]
    num_cols <- sapply(data, is.numeric)
    data[num_cols] <- lapply(data[num_cols], function(x) (x - min(x)) / (max(x) - min(x)))
    # Map categorical columns to integers with even differences
    cat_cols <- sapply(data, function(x) !is.numeric(x) & is.factor(x))
    data[cat_cols] <- lapply(data[cat_cols], function(x) {
        levels_count <- length(levels(x))
        mapped_values <- seq(0, 1, length.out = levels_count)
        mapped_values[match(x, levels(x))]
    })
    res <- stack(data)
    colnames(res) <- c("var_value", "variable")
    return(res[,1])
}


