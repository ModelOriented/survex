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
        sv <- x$result[, cols]
        times <- x$eval_times
        transposed <- as.data.frame(cbind(times = times, sv))
        rownames(transposed) <- NULL
        long_df <- cbind(
            times = transposed$times,
            stack(transposed, select = -times),
            label = label
        )
    })

    transformed_rug_dfs <- lapply(dfl, function(x) {
        label <- attr(x, "label")
        rug_df <- data.frame(times = x$event_times, statuses = as.character(x$event_statuses), label = label)
    })

    rug_df <- do.call(rbind, transformed_rug_dfs)

    long_df <- do.call(rbind, long_df)
    labels <- unique(long_df$label)

    if (!is.null(subtitle) && subtitle == "default") {
        endword <- ifelse(length(labels) > 1, " models", " model")
        subtitle <- paste0("created for the ", paste0(labels, collapse = ", "), endword)
    }

    n_colors <- length(unique(long_df$ind))

    base_plot <- with(long_df, {
        ggplot(data = long_df, aes(x = times, y = values, color = ind)) +
            geom_line(linewidth = 0.8) +
            labs(x = "time", y = "SurvSHAP(t) value", title = title, subtitle = subtitle) +
            xlim(c(0, NA)) +
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
#' @param geom character, one of `"importance"`, `"beeswarm"`, `"profile"` or `"curves"`. Type of chart to be plotted; `"importance"` shows the importance of variables over time and aggregated, `"beeswarm"` shows the distribution of SurvSHAP(t) values for variables and observations, `"profile"` shows the dependence of SurvSHAP(t) values on variable values, `"curves"` shows all SurvSHAP(t) curves for selected variable colored by its value or with functional boxplot if `boxplot = TRUE`.
#' @param ... additional parameters passed to internal functions
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for the XXX model (n = YYY)", where XXX is the explainer label and YYY is the number of observations used for calculations
#' @param max_vars maximum number of variables to be plotted (least important variables are ignored), by default 7
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#'
#' @return An object of the class `ggplot`.
#'
#' @section Plot options:
#'
#' ## `plot.aggregated_surv_shap(geom = "importance")`
#'
#' * `rug` - character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`.
#' * `rug_colors` - character vector containing two colors (containing either hex codes "#FF69B4", or names "blue"). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
#' * `xlab_left, ylab_right` - axis labels for left and right plots (due to different aggregation possibilities)
#'
#'
#' ## `plot.aggregated_surv_shap(geom = "beeswarm")`
#'
#' * no additional parameters
#'
#'
#' ## `plot.aggregated_surv_shap(geom = "profile")`
#'
#' * `variable` - variable for which the profile is to be plotted, by default first from result data
#' * `color_variable` - variable used to denote the color, by default equal to `variable`
#'
#'
#'#' ## `plot.aggregated_surv_shap(geom = "curves")`
#'
#' * `variable` - variable for which SurvSHAP(t) curves are to be plotted, by default first from result data
#' * `boxplot` - whether to plot functional boxplot with marked outliers or all curves colored by variable value
#'
#'
#' @examples
#' \donttest{
#' veteran <- survival::veteran
#' rsf_ranger <- ranger::ranger(
#'     survival::Surv(time, status) ~ .,
#'     data = veteran,
#'     respect.unordered.factors = TRUE,
#'     num.trees = 100,
#'     mtry = 3,
#'     max.depth = 5
#' )
#' rsf_ranger_exp <- explain(
#'     rsf_ranger,
#'     data = veteran[, -c(3, 4)],
#'     y = survival::Surv(veteran$time, veteran$status),
#'     verbose = FALSE
#' )
#'
#' ranger_global_survshap <- model_survshap(
#'     explainer = rsf_ranger_exp,
#'     new_observation = veteran[
#'         c(1:4, 17:20, 110:113, 126:129),
#'         !colnames(veteran) %in% c("time", "status")
#'     ],
#'     y_true = survival::Surv(
#'         veteran$time[c(1:4, 17:20, 110:113, 126:129)],
#'         veteran$status[c(1:4, 17:20, 110:113, 126:129)]
#'     ),
#'     aggregation_method = "integral",
#'     calculation_method = "kernelshap",
#' )
#' plot(ranger_global_survshap)
#' plot(ranger_global_survshap, geom = "beeswarm")
#' plot(ranger_global_survshap, geom = "profile",
#'      variable = "age", color_variable = "karno")
#' plot(ranger_global_survshap, geom = "curves",
#'      variable = "age")
#' plot(ranger_global_survshap, geom = "curves",
#'      variable = "age", boxplot = TRUE)
#' }
#'
#' @export
plot.aggregated_surv_shap <- function(x,
                                      geom = "importance",
                                      ...,
                                      title = "default",
                                      subtitle = "default",
                                      max_vars = 7,
                                      colors = NULL) {

    if (geom == "swarm") {
        geom <- "beeswarm"
    }

    switch(geom,
        "importance" = plot_shap_global_importance(
            x = x,
            ... = ...,
            title = title,
            subtitle = subtitle,
            max_vars = max_vars,
            colors = colors
        ),
        "beeswarm" = plot_shap_global_beeswarm(
            x = x,
            ... = ...,
            title = title,
            subtitle = subtitle,
            max_vars = max_vars,
            colors = colors
        ),
        "profile" = plot_shap_global_profile(
            x = x,
            ... = ...,
            title = title,
            subtitle = subtitle,
            colors = colors
        ),
        "curves" = plot_shap_global_curves(
            x = x,
            ... = ...,
            title = title,
            subtitle = subtitle,
            colors = colors
        ),
        stop("`geom` must be one of 'importance', 'beeswarm', 'profile' or 'curves'")
    )
}

plot_shap_global_importance <- function(x,
                                        ...,
                                        title = "default",
                                        subtitle = "default",
                                        max_vars = 7,
                                        colors = NULL,
                                        rug = "all",
                                        rug_colors = c("#dd0000", "#222222"),
                                        xlab_left = "Average |aggregated SurvSHAP(t)| value",
                                        ylab_right = "Average |SurvSHAP(t)| value") {
    x$result <- aggregate_shap_multiple_observations(x$result, colnames(x$result[[1]]), function(x) mean(abs(x)))
    x$aggregate <- apply(do.call(rbind, x$aggregate), 2, function(x) mean(abs(x)))

    right_plot <- plot.surv_shap(
        x = x,
        title = NULL,
        subtitle = NULL,
        max_vars = max_vars,
        colors = NULL,
        rug = rug,
        rug_colors = rug_colors
    ) +
        labs(y = ylab_right) +
        theme_default_survex()

    label <- attr(x, "label")
    long_df <- stack(x$aggregate)
    long_df <- long_df[order(long_df$values, decreasing = TRUE), ][1:min(max_vars, length(x$aggregate)), ]

    if (!is.null(title) && title == "default") {
        title <- "Feature importance according to aggregated |SurvSHAP(t)|"
    }
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0(
            "created for the ", label, " model ",
            "(n=", x$n_observations, ")"
        )
    }

    if (is.null(colors)) {
        colors <- c(
            low = "#9fe5bd",
            mid = "#46bac2",
            high = "#371ea3"
        )
    }

    left_plot <- with(long_df, {
        ggplot(long_df, aes(x = values, y = reorder(ind, values))) +
            geom_col(fill = colors[2]) +
            theme_default_survex() +
            labs(x = xlab_left, y = "variable") +
            theme(axis.title.y = element_blank())
    })


    pl <- left_plot +
        right_plot +
        patchwork::plot_layout(widths = c(3, 5), guides = "collect") +
        patchwork::plot_annotation(title = title, subtitle = subtitle) &
        theme_default_survex() &
        theme(
            legend.position = "top",
            plot.title = element_text(color = "#371ea3", size = 16, hjust = 0),
            plot.subtitle = element_text(color = "#371ea3", hjust = 0),
        )

    return(pl)
}

plot_shap_global_beeswarm <- function(x,
                                      ...,
                                      title = "default",
                                      subtitle = "default",
                                      max_vars = 7,
                                      colors = NULL) {
    df <- as.data.frame(do.call(rbind, x$aggregate))
    cols <- names(sort(colMeans(abs(df))))[1:min(max_vars, length(df))]
    df <- df[, cols]
    df <- stack(df)
    colnames(df) <- c("shap_value", "variable")

    original_values <- as.data.frame(x$variable_values)[, cols]
    var_value <- preprocess_values_to_common_scale(original_values)
    df <- cbind(df, var_value)

    label <- attr(x, "label")
    if (!is.null(title) && title == "default") {
        title <- "Aggregated SurvSHAP(t) values summary"
    }
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0(
            "created for the ", label, " model ",
            "(n=", x$n_observations, ")"
        )
    }

    if (is.null(colors)) {
        colors <- c(
            low = "#9fe5bd",
            mid = "#46bac2",
            high = "#371ea3"
        )
    }

    with(df, {
        ggplot(data = df, aes(x = shap_value, y = variable, color = var_value)) +
            geom_vline(xintercept = 0, color = "#ceced9", linetype = "solid") +
            geom_jitter(width = 0, height = 0.15) +
            scale_color_gradient2(
                name = "Variable value",
                low = colors[1],
                mid = colors[2],
                high = colors[3],
                midpoint = 0.5,
                limits = c(0, 1),
                breaks = c(0, 1),
                labels = c("", "")
            ) +
            labs(
                title = title, subtitle = subtitle,
                x = "Aggregated SurvSHAP(t) value",
                y = "Variable"
            ) +
            theme_default_survex() +
            theme(legend.position = "bottom") +
            guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))
    })
}

plot_shap_global_profile <- function(x,
                                     ...,
                                     variable = NULL,
                                     color_variable = NULL,
                                     title = "default",
                                     subtitle = "default",
                                     colors = NULL) {
    df <- as.data.frame(do.call(rbind, x$aggregate))

    if (is.null(variable)) {
        variable <- colnames(df)[1]
        warning("`variable` was not specified, the first from the result will be used.")
    }
    if (is.null(color_variable)) {
        color_variable <- variable
        warning("`color_variable` was not specified, the first from the result will be used.")
    }

    shap_val <- df[, variable]

    original_values <- as.data.frame(x$variable_values)
    var_vals <- original_values[, c(variable, color_variable)]

    df <- cbind(shap_val, var_vals)
    colnames(df) <- c("shap_val", "variable_val", "color_variable_val")

    label <- attr(x, "label")
    if (!is.null(title) && title == "default") {
        title <- "Aggregated SurvSHAP(t) profile"
    }
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", unique(variable), " variable")
    }

    p <- with(df, {
        ggplot(df, aes(x = variable_val, y = shap_val, color = color_variable_val)) +
            geom_hline(yintercept = 0, color = "#ceced9", linetype = "solid") +
            geom_point() +
            geom_rug(aes(x = variable_val), inherit.aes = F, color = "#ceced9") +
            labs(
                x = paste(variable, "value"),
                y = "Aggregated SurvSHAP(t) value",
                title = title,
                subtitle = subtitle
            ) +
            theme_default_survex() +
            theme(legend.position = "bottom")
    })

    if (is.factor(df$color_variable_val) || is.character(df$color_variable_val)) {
        p + scale_color_manual(
            name = paste(color_variable, "value"),
            values = generate_discrete_color_scale(length(unique(df$color_variable_val)), colors)
        )
    } else {
        if (is.null(colors)) {
            colors <- c(
                low = "#9fe5bd",
                mid = "#46bac2",
                high = "#371ea3"
            )
        }
        p + scale_color_gradient2(
            name = paste(color_variable, "value"),
            low = colors[1],
            mid = colors[2],
            high = colors[3],
            midpoint = median(df$color_variable_val)
        )
    }
}


plot_shap_global_curves <- function(x,
                                    ...,
                                    variable = NULL,
                                    boxplot = FALSE,
                                    coef = 1.5,
                                    title = "default",
                                    subtitle = "default",
                                    colors = NULL,
                                    rug = "all",
                                    rug_colors = c("#dd0000", "#222222")) {

    if (is.null(variable)) {
        variable <- colnames(df)[1]
        warning("`variable` was not specified, the first from the result will be used.")
    }

    label <- attr(x, "label")
    if (!is.null(title) && title == "default") {
        title <- "SurvSHAP(t) curves"
    }
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", unique(variable), " variable")
    }

    if (!boxplot){
        df <- as.data.frame(do.call(rbind, x$result))
        df$obs <- rep(1:x$n_observations,
                      each = length(x$eval_times))
        df$time <- rep(x$eval_times,
                       times = x$n_observations)
        df$varval <- rep(x$variable_values[[variable]],
                         each = length(x$eval_times))

        if (is.null(colors) || length(colors) < 3) {
            colors <- c(
                low = "#9fe5bd",
                mid = "#46bac2",
                high = "#371ea3"
            )
        }
        base_plot <- with(df,
                          {ggplot(df, aes(x = time, y = !!sym(variable), color = varval)) +
                            geom_hline(yintercept = 0, alpha = 0.5, color = "black") +
                            geom_line(aes(group = obs), alpha = 0.5) +
                            theme_default_survex() +
                            theme(legend.position = "bottom") +
                            labs(y = "SurvSHAP(t) value",
                                 title = title,
                                 subtitle = subtitle)
        })

        if (is.factor(x$variable_values[[variable]]) || is.character(x$variable_values[[variable]])) {
            base_plot <- base_plot + scale_color_manual(
                name = paste(variable, "value"),
                values = generate_discrete_color_scale(length(unique(df$varval)), colors)
            )
        } else {
            base_plot <- base_plot + scale_color_gradient2(
                name = paste(variable, "value"),
                low = colors[1],
                mid = colors[2],
                high = colors[3],
                midpoint = median(as.numeric(as.character(df$varval)))
            )
        }
    } else {
        if (is.null(colors) || length(colors) < 3) {
            colors <- c(
                background = "#000000",
                boxplot = "#9fe5bd",
                outliers = "#371ea3"
            )
        }

        df <- t(sapply( x$result, function(x) x[[variable]]))
        n <- x$n_observations
        p <- dim(df)[2]

        rmat <- apply(df, 2, rank)
        down <- rmat-1
        up <- n-rmat
        depth <- (rowSums(up*down)/p+n-1)/choose(n,2)

        index <- order(depth,decreasing=TRUE)
        median <- colMeans(df[which(depth==max(depth)), ,drop=FALSE])

        m <- ceiling(n * 0.5)
        central_region <- df[index[1:m], ]
        outer_region <- df[index[(m+1):n], ]
        lower_bound <- apply(central_region, 2, min)
        upper_bound <- apply(central_region, 2, max)

        iqr <- upper_bound - lower_bound
        outlier_indices <- which(colSums((t(df) <= median - coef * iqr) +
                                        (t(df) >= median + coef * iqr)) > 0)

        if (length(outlier_indices) > 0){
            nonoutliers <- df[-outlier_indices,]
            outliers <- df[outlier_indices, ]
            df_outliers <- data.frame(x = x$eval_times,
                                      y = as.vector(t(outliers)),
                                      obs = rep(1:length(outliers), each = length(x$eval_times)))
        } else {
            nonoutliers <- df
            df_outliers <- NULL
        }

        whisker_lower_bound <- apply(nonoutliers, 2, min)
        whisker_upper_bound <- apply(nonoutliers, 2, max)

        df_all <- data.frame(x = x$eval_times,
                             y = as.vector(t(df)),
                             obs = rep(1:n, each = length(x$eval_times)))

        df_res <- data.frame(x = x$eval_times,
                             median = median,
                             lower_bound = lower_bound,
                             upper_bound = upper_bound,
                             whisker_lower_bound = whisker_lower_bound,
                             whisker_upper_bound = whisker_upper_bound)

        base_plot <- with(list(df_all, df_outliers, df_res),
             {ggplot() +
                 geom_hline(yintercept = 0, alpha = 0.5, color = colors[1]) +
                 geom_line(data = df_all, aes(x = x, y = y, group = obs), col = colors[1], alpha = 0.2, linewidth = 0.2) +
                 geom_ribbon(data = df_res, aes(x = x, ymin = lower_bound, ymax = upper_bound),
                             fill = colors[1], alpha = 0.4) +
                 geom_line(data = df_res, aes(x = x, y = median), col = colors[2], linewidth = 1) +
                 geom_line(data = df_res, aes(x = x, y = lower_bound), col = colors[2], linewidth = 1) +
                 geom_line(data = df_res, aes(x = x, y = upper_bound), col = colors[2], linewidth = 1) +
                 geom_line(data = df_res, aes(x = x, y = whisker_lower_bound), col = colors[2], lty = 2, linewidth = 1) +
                 geom_line(data = df_res, aes(x = x, y = whisker_upper_bound), col = colors[2], lty = 2, linewidth = 1) +
                 theme_default_survex() +
                 labs(x = "time",
                      y = "SurvSHAP(t) value",
                      title = title,
                      subtitle = subtitle)
            })
        if (length(outlier_indices) > 0){
            base_plot <- base_plot +
                geom_line(data = df_outliers, aes(x = x, y = y, group = obs),
                          col = colors[3], linewidth = 0.5, alpha = 0.1)
            cat("Observations with outlying SurvSHAP(t) values:\n")
            print(x$variable_values[outlier_indices,])
        }
    }

    rug_df <- data.frame(times = x$event_times, statuses = as.character(x$event_statuses))
    return_plot <- add_rug_to_plot(base_plot, rug_df, rug, rug_colors)
    return(return_plot)
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
    return(res[, 1])
}
