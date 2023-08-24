#' Plot Predict Profile for Survival Models
#'
#' This function plots objects of class `"predict_profile_survival"` created using
#' the `predict_profile()` function.
#'
#' @param x an object of class `predict_profile_survival` to be plotted
#' @param ... additional objects of class `"predict_profile_survival"` to be plotted together. Only available for `geom = "time"`.
#' @param geom character, either `"time"` or `"variable"`. Selects the type of plot to be prepared. If `"time"` then the x-axis represents survival times, and variable is denoted by colors, if `"variable"` then the x-axis represents the variable values, and y-axis represents the predictions at selected time points.
#' @param variables character, names of the variables to be plotted. When `geom = "variable"` it needs to be a name of a single variable, when `geom = "time"` it can be a vector of variable names. If `NULL` (default) then first variable (for `geom = "variable"`) or all variables (for `geom = "time"`) are plotted.
#' @param variable_type character, either `"numerical"`, `"categorical"` or `NULL` (default), select only one type of variable for plotting, or leave `NULL` for all.  Only used when `geom = "time"`.
#' @param facet_ncol number of columns for arranging subplots.  Only used when `geom = "time"`.
#' @param numerical_plot_type character, either `"lines"`, or `"contours"` selects the type of numerical variable plots.  Only used when `geom = "time"`.
#' @param times numeric vector, times for which the profile should be plotted, the times must be present in the 'times' field of the explainer. If `NULL` (default) then the median survival time (if available) or the median time from the explainer object is used. Only used when `geom = "variable"` and `marginalize_over_time = FALSE`.
#' @param marginalize_over_time logical, if `TRUE` then the profile is calculated for all times and then averaged over time, if `FALSE` (default) then the profile is calculated for each time separately. Only used when `geom = "variable"`.
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' @param rug character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`.
#' @param rug_colors character vector containing two colors (containing either hex codes `"#FF69B4"`, or names `"blue"`). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
#'
#' @return A collection of `ggplot` objects arranged with the `patchwork` package.
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
#'
#' @export
plot.predict_profile_survival <- function(x,
                                          ...,
                                          geom = "time",
                                          variables = NULL,
                                          variable_type = NULL,
                                          facet_ncol = NULL,
                                          numerical_plot_type = "lines",
                                          times = NULL,
                                          marginalize_over_time = FALSE,
                                          title = "default",
                                          subtitle = "default",
                                          colors = NULL,
                                          rug = "all",
                                          rug_colors = c("#dd0000", "#222222")) {

    if (!geom %in% c("time", "variable")) {
        stop("`geom` needs to be one of 'time' or 'variable'.")
    }

    if (!(numerical_plot_type %in% c("lines", "contours"))) {
        stop("`numerical_plot_type` needs to be 'lines' or 'contours'")
    }

    if (!is.null(variable_type) &&
        !(variable_type %in% c("numerical", "categorical"))) {
        stop("`variable_type` needs to be 'numerical' or 'categorical'")
    }

    if (title == "default"){
        title <- "Ceteris paribus survival profile"
    }

    if (geom == "variable") {
        pl <- plot2_cp(
            x = x,
            variable = variables,
            times = times,
            marginalize_over_time = marginalize_over_time,
            ... = ...,
            title = title,
            subtitle = subtitle,
            colors = colors
        )
        if (x$center) {
            pl <- pl + labs(y = "centered profile value")
        }
        return(pl)
    }

    lapply(list(x, ...), function(x) {
        if (!inherits(x, "predict_profile_survival")) {
            stop("All ... must be objects of class `predict_profile_survival`.")
        }
    })
    explanations_list <- c(list(x), list(...))
    num_models <- length(explanations_list)

    if (num_models == 1) {
        result <- prepare_ceteris_paribus_plots(
            x,
            colors,
            variable_type,
            facet_ncol,
            variables,
            numerical_plot_type,
            title,
            subtitle,
            rug,
            rug_colors
        )
        return(result)
    }

    return_list <- list()
    labels <- list()
    for (i in 1:num_models) {
        this_title <- unique(explanations_list[[i]]$result$`_label_`)
        return_list[[i]] <- prepare_ceteris_paribus_plots(
            explanations_list[[i]],
            colors,
            variable_type,
            1,
            variables,
            numerical_plot_type,
            this_title,
            NULL,
            rug,
            rug_colors
        )
        labels[[i]] <- c(this_title, rep("", length(return_list[[i]]$patches) - 2))
    }

    labels <- unlist(labels)
    return_plot <- patchwork::wrap_plots(return_list, nrow = 1, tag_level = "keep") +
        patchwork::plot_annotation(title, tag_levels = list(labels)) & theme_default_survex()

    return(return_plot)
}

#' @keywords internal
plot2_cp <- function(x,
                     variable,
                     times = NULL,
                     marginalize_over_time = FALSE,
                     ...,
                     title = "default",
                     subtitle = "default",
                     colors = NULL) {
    if (is.null(variable)) {
        variable <- unique(x$result$`_vname_`)[1]
        warning("Plot will be prepared for the first variable from the explanation `result`. \nFor another variable, set the value of `variable`.")
    }

    if (!is.character(variable)) {
        stop("The variable must be specified by name")
    }

    if (length(variable) > 1) {
        stop("Only one variable can be specified for `geom`='variable'")
    }

    if (!variable %in% x$result$`_vname_`) {
        stop(paste0("Variable ", variable, " not found"))
    }

    if (is.null(times)) {
        if (marginalize_over_time){
            times <- x$eval_times
            warning("Plot will be prepared with marginalization over all time points from the explainer's `times` vector. \nFor subset of time points, set the value of `times`.")
        } else if (!is.null(x$median_survival_time)){
            times <- x$median_survival_time
            warning("Plot will be prepared for the median survial time. For another time point, set the value of `times`.")
        } else {
            times <- quantile(x$eval_times, p = 0.5, type = 1)
            warning("Plot will be prepared for the median time point from the explainer's `times` vector. For another time point, set the value of `times`.")
        }
    }

    if (!all(times %in% x$eval_times)) {
        stop(paste0(
            "For one of the provided times the explanations has not been calculated not found.
         Please modify the times argument in your explainer or use only values from the following:  ",
            paste(x$eval_times, collapse = ", ")
        ))
    }

    single_timepoint <- ((length(times) == 1) || marginalize_over_time)
    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", unique(variable), " variable")
        if (single_timepoint && !marginalize_over_time) {
            subtitle <- paste0(subtitle, " and time =", times)
        }
    }

    is_categorical <- (unique(x$result[x$result$`_vname_` == variable, "_vtype_"]) == "categorical")

    if (single_timepoint) {
        ice_df <- x$result[(x$result$`_vname_` == variable) & (x$result$`_times_` %in% times), c(variable, "_yhat_")]
        colnames(ice_df) <- c(variable, "ice")
    } else {
        ice_df <- x$result[(x$result$`_vname_` == variable) & (x$result$`_times_` %in% times), c(variable, "_times_", "_yhat_")]
        colnames(ice_df) <- c(variable, "time", "ice")
        ice_df$time <- as.factor(ice_df$time)
    }

    if (is_categorical) {
        ice_df[, variable] <- as.factor(ice_df[, variable])
    }

    feature_name_sym <- sym(variable)

    if (marginalize_over_time) {
        ice_df <- aggregate(ice ~ ., data = ice_df, mean)
        color_scale <- generate_discrete_color_scale(1, colors)
    } else {
        if (is.null(colors) || length(colors) < 3) {
            color_scale <- c(
                low = "#9fe5bd",
                mid = "#46bac2",
                high = "#371ea3"
            )
        }
    }

    if (is_categorical) {
        pl <- plot_ice_cat(
            ice_df = ice_df,
            feature_name_sym = feature_name_sym,
            single_timepoint = single_timepoint,
            colors = color_scale
        )
    } else {
        ice_df[, 1] <- as.numeric(as.character(ice_df[, 1]))
        x_width <- diff(range(ice_df[, variable]))
        pl <- plot_ice_num(
            ice_df = ice_df,
            feature_name_sym = feature_name_sym,
            x_width = x_width,
            single_timepoint = single_timepoint,
            colors = color_scale
        )
    }

    pl +
        labs(
            title = title,
            subtitle = subtitle
        ) +
        theme_default_survex()
}



plot_ice_num <-function(ice_df,
                        feature_name_sym,
                        x_width = x_width,
                        single_timepoint,
                        colors){
    with(ice_df, {
        if (single_timepoint == TRUE) {
                ggplot(data = ice_df, aes(x = !!feature_name_sym, y = ice)) +
                    geom_line(color = colors_discrete_drwhy(1), linewidth=1)
        } else {
            ice_df$time <- as.numeric(as.character(ice_df$time))
            ggplot(data = ice_df, aes(x = !!feature_name_sym, y = ice)) +
                geom_line(aes(color = time, group = time), linewidth=1) +
                scale_colour_gradient2(
                    low = colors[1],
                    mid = colors[2],
                    high = colors[3],
                    midpoint = median(as.numeric(as.character(pdp_dt$time)))
                )
        }
    })
}


plot_ice_cat <- function(ice_df,
                         feature_name_sym,
                         single_timepoint,
                         colors){
    with(ice_df, {
        if (single_timepoint == TRUE) {
            ggplot(data = ice_df, aes(x = !!feature_name_sym, y = ice), ) +
                geom_bar(stat = "identity", width = 0.5, fill = colors_discrete_drwhy(1)) +
                scale_y_continuous() +
                scale_fill_manual(name = "time", values = colors) +
                geom_hline(yintercept = 0, linetype="dashed")
        } else {
            ggplot(data = ice_df, aes(x = !!feature_name_sym, y = ice, fill = time)) +
                geom_bar(stat = "identity", width = 0.5, position = "dodge") +
                scale_y_continuous() +
                scale_fill_manual(name = "time", values = colors)
        }
    })
}
