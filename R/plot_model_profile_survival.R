#' Plot Model Profile for Survival Models
#'
#' This function plots objects of class `"model_profile_survival"` created
#' using the `model_profile()` function.
#'
#' @param x an object of class `model_profile_survival` to be plotted
#' @param ... additional objects of class `model_profile_survival` to be plotted together. Only available for `geom = "time"`.
#' @param geom character, either `"time"` or `"variable"`. Selects the type of plot to be prepared. If `"time"` then the x-axis represents survival times, and variable is denoted by colors, if `"variable"` then the x-axis represents the variable values, and y-axis represents the predictions at selected time points.
#' @param variables character, names of the variables to be plotted. When `geom = "variable"` it needs to be a name of a single variable, when `geom = "time"` it can be a vector of variable names. If `NULL` (default) then first variable (for `geom = "variable"`) or all variables (for `geom = "time"`) are plotted.
#' @param variable_type character, either `"numerical"`, `"categorical"` or `NULL` (default), select only one type of variable for plotting, or leave `NULL` for all. Only used when `geom = "time"`.
#' @param facet_ncol number of columns for arranging subplots. Only used when `geom = "time"`.
#' @param numerical_plot_type character, either `"lines"`, or `"contours"` selects the type of numerical variable plots. Only used when `geom = "time"`.
#' @param times numeric vector, times for which the profile should be plotted, the times must be present in the 'times' field of the explainer. If `NULL` (default) then the median survival time (if available) or the median time from the explainer object is used. Only used when `geom = "variable"` and `marginalize_over_time = FALSE`.
#' @param marginalize_over_time logical, if `TRUE` then the profile is calculated for all times and then averaged over time, if `FALSE` (default) then the profile is calculated for each time separately. Only used when `geom = "variable"`.
#' @param plot_type character, one of `"pdp"`, `"ice"`, `"pdp+ice"`, or `NULL` (default). If `NULL` then the type of plot is chosen automatically based on the number of variables to be plotted. Only used when `geom = "variable"`.
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `"default"` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue").
#' @param rug character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`. Only used when `geom = "time"`.
#' @param rug_colors character vector containing two colors (containing either hex codes "#FF69B4", or names "blue"). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
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
#' m_prof <- model_profile(exp, categorical_variables = "trt")
#'
#' plot(m_prof)
#'
#' plot(m_prof, numerical_plot_type = "contours")
#'
#' plot(m_prof, variables = c("trt", "age"), facet_ncol = 1)
#'
#' plot(m_prof, geom = "variable", variables = "karno", plot_type = "pdp+ice")
#'
#' plot(m_prof, geom = "variable", times = exp$times[c(5, 10)],
#'      variables = "karno", plot_type = "pdp+ice")
#'
#' plot(m_prof, geom = "variable", times = exp$times[c(5, 10)],
#'      variables = "trt", plot_type = "pdp+ice")
#' }
#'
#' @export
plot.model_profile_survival <- function(x,
                                        ...,
                                        geom = "time",
                                        variables = NULL,
                                        variable_type = NULL,
                                        facet_ncol = NULL,
                                        numerical_plot_type = "lines",
                                        times = NULL,
                                        marginalize_over_time = FALSE,
                                        plot_type = NULL,
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

    if (!is.null(title) && title == "default") {
        if (x$type == "partial") {
            title <- "Partial dependence survival profiles"
            if (geom == "variable") {
                title <- "default"
            }
        }
        if (x$type == "accumulated") {
            title <- "Accumulated local effects survival profiles"
        }
    }

    if (geom == "variable") {
        pl <- plot2_mp(
            x = x,
            variable = variables,
            times = times,
            marginalize_over_time = marginalize_over_time,
            plot_type = plot_type,
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
        if (!inherits(x, "model_profile_survival")) {
            stop("All ... must be objects of class `model_profile_survival`.")
        }
    })
    explanations_list <- c(list(x), list(...))
    num_models <- length(explanations_list)

    if (num_models == 1) {
        result <- prepare_model_profile_plots(
            x,
            variables = variables,
            variable_type = variable_type,
            facet_ncol = facet_ncol,
            numerical_plot_type = numerical_plot_type,
            title = title,
            subtitle = subtitle,
            colors = colors,
            rug = rug,
            rug_colors = rug_colors
        )
        return(result)
    }

    return_list <- list()
    labels <- list()
    for (i in 1:num_models) {
        this_title <- unique(explanations_list[[i]]$result$`_label_`)
        return_list[[i]] <- prepare_model_profile_plots(
            explanations_list[[i]],
            variables = variables,
            variable_type = variable_type,
            facet_ncol = 1,
            numerical_plot_type = numerical_plot_type,
            title = this_title,
            subtitle = NULL,
            colors = colors,
            rug = rug,
            rug_colors = rug_colors
        )
        labels[[i]] <- c(this_title, rep("", length(return_list[[i]]$patches) - 2))
    }

    labels <- unlist(labels)
    return_plot <- patchwork::wrap_plots(return_list, nrow = 1, tag_level = "keep") +
        patchwork::plot_annotation(title, tag_levels = list(labels)) & theme_default_survex()

    return(return_plot)
}


#' @keywords internal
plot2_mp <- function(x,
                  variable,
                  times = NULL,
                  marginalize_over_time = FALSE,
                  plot_type = NULL,
                  ...,
                  title = "default",
                  subtitle = "default",
                  colors = NULL) {
    if (is.null(plot_type)) {
        if (x$type == "accumulated") {
            plot_type <- "ale"
        } else if (x$type == "partial") plot_type <- "pdp+ice"
    }

    if (plot_type == "ice") {
        title <- "Individual conditional expectation survival profiles"
    } else if (plot_type == "pdp+ice") {
        title <- "Partial dependence with individual conditional expectation survival profiles"
    } else if (plot_type == "pdp") {
        title <- "Partial dependence survival profiles"
    }

    if (x$type == "accumulated" && plot_type != "ale") {
        stop("For accumulated local effects explanations only plot_type = 'ale' is available")
    }

    if (!plot_type %in% c("pdp", "ice", "pdp+ice", "ale")) {
        stop("plot_type must be one of 'pdp', 'ice', 'pdp+ice'")
    }

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
    ice_needed <- plot_type %in% c("pdp+ice", "ice")

    if (single_timepoint) {
        pdp_df <- x$result[(x$result$`_vname_` == variable) & (x$result$`_times_` %in% times), c("_x_", "_yhat_")]
        colnames(pdp_df) <- c(variable, "pd")
    } else {
        pdp_df <- x$result[(x$result$`_vname_` == variable) & (x$result$`_times_` %in% times), c("_x_", "_times_", "_yhat_")]
        colnames(pdp_df) <- c(variable, "time", "pd")
        pdp_df$time <- as.factor(pdp_df$time)
    }

    ice_df <- NULL
    if (ice_needed) {
        ice_df <- x$cp_profiles$result[(x$cp_profiles$result$`_vname_` == variable) &
            (x$cp_profiles$result$`_times_` %in% times), ]

        if (single_timepoint) {
            ice_df$`_times_` <- NULL
        } else {
            colnames(ice_df)[colnames(ice_df) == "_times_"] <- "time"
            ice_df$time <- as.factor(ice_df$time)
        }

        if (is_categorical) {
            ice_df[, variable] <- as.factor(ice_df[, variable])
        }

        ice_df$`_vname_` <- NULL
        ice_df$`_vtype_` <- NULL
        ice_df$`_label_` <- NULL

        colnames(ice_df)[colnames(ice_df) == "_ids_"] <- "id"
        colnames(ice_df)[colnames(ice_df) == "_yhat_"] <- "predictions"

        y_floor_ice <- floor(min(ice_df[, "predictions"]) * 10) / 10
        y_ceiling_ice <- ceiling(max(ice_df[, "predictions"]) * 10) / 10
    }


    if (is_categorical) {
        pdp_df[, variable] <- as.factor(pdp_df[, variable])
    }



    feature_name_sym <- sym(variable)

    data_df <- x$cp_profiles$variable_values

    y_floor_pd <- floor(min(pdp_df[, "pd"]) * 10) / 10
    y_ceiling_pd <- ceiling(max(pdp_df[, "pd"]) * 10) / 10

    if (marginalize_over_time) {
        pdp_df <- aggregate(pd ~ ., data = pdp_df, mean)
        ice_df <- aggregate(predictions ~ ., data = ice_df, mean)
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
        pl <- plot_pdp_cat(
            pdp_dt = pdp_df,
            ice_dt = ice_df,
            data_dt = data_df,
            feature_name_count_sym = feature_name_sym,
            y_floor_ice = y_floor_ice,
            y_ceiling_ice = y_ceiling_ice,
            y_floor_pd = y_floor_pd,
            y_ceiling_pd = y_ceiling_pd,
            plot_type = plot_type,
            single_timepoint = single_timepoint,
            colors = color_scale
        )
    } else {
        pdp_df[, 1] <- as.numeric(as.character(pdp_df[, 1]))
        x_width <- diff(range(pdp_df[, variable]))
        pl <- plot_pdp_num(
            pdp_dt = pdp_df,
            ice_dt = ice_df,
            data_dt = data_df,
            feature_name_sym = feature_name_sym,
            y_floor_ice = y_floor_ice,
            y_ceiling_ice = y_ceiling_ice,
            y_floor_pd = y_floor_pd,
            y_ceiling_pd = y_ceiling_pd,
            x_width = x_width,
            plot_type = plot_type,
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

plot_pdp_num <- function(pdp_dt,
                         ice_dt,
                         data_dt,
                         feature_name_sym,
                         y_floor_ice,
                         y_ceiling_ice,
                         y_floor_pd,
                         y_ceiling_pd,
                         x_width,
                         plot_type,
                         single_timepoint,
                         colors) {
    with(pdp_dt, { # to get rid of Note: no visible binding for global variable ...
        if (single_timepoint == TRUE) { ## single timepoint
            if (plot_type == "ice") {
                ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                    geom_line(alpha = 0.2, mapping = aes(group = id), color = colors_discrete_drwhy(1)) +
                    geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = position_jitter(width = 0.01 * x_width)) +
                    ylim(y_floor_ice, y_ceiling_ice)
            }
            # PDP + ICE
            else if (plot_type == "pdp+ice") {
                ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                    geom_line(mapping = aes(group = id), alpha = 0.2) +
                    geom_line(data = pdp_dt, aes(x = !!feature_name_sym, y = pd), linewidth = 2, color = colors_discrete_drwhy(1)) +
                    geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = position_jitter(width = 0.01 * x_width)) +
                    ylim(y_floor_ice, y_ceiling_ice)
            }
            # PDP
            else if (plot_type == "pdp" || plot_type == "ale") {
                ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
                    geom_line(color = colors_discrete_drwhy(1)) +
                    geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = position_jitter(width = 0.01 * x_width)) +
                    ylim(y_floor_pd, y_ceiling_pd)
            }
        } else { ## multiple time points
            pdp_dt$time <- as.numeric(as.character(pdp_dt$time))
            if (!is.null(ice_dt)) {
                ice_dt$time <- as.numeric(as.character(ice_dt$time))
            }

            if (plot_type == "ice") {
                ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                    geom_line(alpha = 0.2, mapping = aes(group = interaction(id, time), color = time)) +
                    geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides = "b", alpha = 0.8, position = position_jitter(width = 0.01 * x_width)) +
                    scale_colour_gradient2(
                        low = colors[1],
                        mid = colors[2],
                        high = colors[3],
                        midpoint = median(as.numeric(as.character(pdp_dt$time)))
                    ) +
                    ylim(y_floor_ice, y_ceiling_ice)
            }
            # PDP + ICE
            else if (plot_type == "pdp+ice") {
                ggplot() +
                    geom_line(data = ice_dt, aes(x = !!feature_name_sym, y = predictions, group = interaction(id, time), color = time), alpha = 0.1) +
                    geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, color = time, group = time), linewidth = 1.5, lineend = "round", linejoin = "round") +
                    geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, group = time), color = "black", linewidth = 0.5, linetype = "dashed", lineend = "round", linejoin = "round") +
                    geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides = "b", alpha = 0.8, position = position_jitter(width = 0.01 * x_width)) +
                    scale_colour_gradient2(
                        low = colors[1],
                        mid = colors[2],
                        high = colors[3],
                        midpoint = median(as.numeric(as.character(pdp_dt$time)))
                    ) +
                    ylim(y_floor_ice, y_ceiling_ice)
            }
            # PDP
            else if (plot_type == "pdp" || plot_type == "ale") {
                ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
                    geom_line(aes(color = time, group = time)) +
                    geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = position_jitter(width = 0.01 * x_width)) +
                    scale_colour_gradient2(
                        low = colors[1],
                        mid = colors[2],
                        high = colors[3],
                        midpoint = median(as.numeric(as.character(pdp_dt$time)))
                    ) +
                    ylim(y_floor_pd, y_ceiling_pd)
            }
        }
    })
}

plot_pdp_cat <- function(pdp_dt,
                         ice_dt,
                         data_dt,
                         feature_name_count_sym,
                         y_floor_ice,
                         y_ceiling_ice,
                         y_floor_pd,
                         y_ceiling_pd,
                         plot_type,
                         single_timepoint,
                         colors) {
    with(pdp_dt, { # to get rid of Note: no visible binding for global variable ...
        if (single_timepoint == TRUE) { ## single timepoint
            if (plot_type == "ice") {
                ggplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions)) +
                    geom_boxplot(alpha = 0.2, color = colors_discrete_drwhy(1)) +
                    scale_color_manual(name = "time", values = colors)
            } else if (plot_type == "pdp+ice") {
                ggplot() +
                    geom_boxplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions), alpha = 0.2, color = colors_discrete_drwhy(1), width = 0.7) +
                    geom_line(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, group = 1), linewidth = 2, color = colors_discrete_drwhy(1), position = position_dodge(0.7)) +
                    scale_color_manual(name = "time", values = colors)
            } else if (plot_type == "pdp" || plot_type == "ale") {
                ggplot(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd), ) +
                    geom_bar(stat = "identity", width = 0.5, fill = colors_discrete_drwhy(1)) +
                    scale_y_continuous(expand = c(0, NA)) +
                    scale_fill_manual(name = "time", values = colors)
            }
        } else {
            pdp_dt$time <- as.numeric(as.character(pdp_dt$time))
            if (!is.null(ice_dt)) {
                ice_dt$time <- as.numeric(as.character(ice_dt$time))
            }
            if (plot_type == "ice") {
                ggplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions)) +
                    geom_boxplot(mapping = aes(color = time, group = interaction(time, !!feature_name_count_sym)), alpha = 0.2) +
                    scale_colour_gradient2(
                        low = colors[1],
                        mid = colors[2],
                        high = colors[3],
                        midpoint = median(as.numeric(as.character(ice_dt$time)))
                    )
            } else if (plot_type == "pdp+ice") {
                ggplot(mapping = aes(color = time)) +
                    geom_boxplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions, group = interaction(time, !!feature_name_count_sym)), alpha = 0.2, width = 0.7) +
                    geom_line(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, group = time), linewidth = 0.6, position = position_dodge(0.7)) +
                    scale_colour_gradient2(
                        low = colors[1],
                        mid = colors[2],
                        high = colors[3],
                        midpoint = median(as.numeric(as.character(pdp_dt$time)))
                    )
            } else if (plot_type == "pdp" || plot_type == "ale") {
                ggplot(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, fill = time, group = time)) +
                    geom_bar(stat = "identity", width = 0.5, position = "dodge") +
                    scale_y_continuous(expand = c(0, NA)) +
                    scale_fill_gradient2(
                        low = colors[1],
                        mid = colors[2],
                        high = colors[3],
                        midpoint = median(as.numeric(as.character(pdp_dt$time)))
                    )
            }
        }
    })
}




prepare_model_profile_plots <- function(x,
                                        variables = NULL,
                                        variable_type = NULL,
                                        facet_ncol = NULL,
                                        numerical_plot_type = "lines",
                                        title = "default",
                                        subtitle = "default",
                                        colors = NULL,
                                        rug = rug,
                                        rug_colors = rug_colors) {
    rug_df <- data.frame(times = x$event_times, statuses = as.character(x$event_statuses), label = unique(x$result$`_label_`))
    aggregated_profiles <- x$result
    class(aggregated_profiles) <- "data.frame"

    all_variables <- unique(aggregated_profiles$`_vname_`)
    if (!is.null(variables)) {
        all_variables <- intersect(all_variables, variables)
        if (length(all_variables) == 0) {
            stop(paste0(
                "variables do not overlap with ",
                paste(all_variables, collapse = ", ")
            ))
        }
    }
    aggregated_profiles <-
        aggregated_profiles[aggregated_profiles$`_vname_` %in% all_variables, ]


    if (!is.null(subtitle) && subtitle == "default") {
        labels <- unique(aggregated_profiles$`_label_`)
        endword <- ifelse(length(labels) > 1, " models", " model")
        subtitle <- paste0("created for the ", paste0(labels, collapse = ", "), endword)
    }

    if (is.null(variables)) {
        variables <- unique(aggregated_profiles$`_vname_`)
    }

    if (!is.null(variable_type) && variable_type == "numerical") {
        aggregated_profiles <- aggregated_profiles[aggregated_profiles$`_vtype_` == "numerical", ]
    }

    if (!is.null(variable_type) && variable_type == "categorical") {
        aggregated_profiles <- aggregated_profiles[aggregated_profiles$`_vtype_` == "categorical", ]
    }

    variables <- intersect(variables, unique(aggregated_profiles$`_vname_`))
    n_colors <- length(unique(aggregated_profiles$`_x_`))

    aggregated_profiles$`_real_point_` <- FALSE

    pl <- plot_individual_ceteris_paribus_survival(aggregated_profiles, variables, colors, numerical_plot_type, rug_df, rug, rug_colors, x$center, x$output_type)

    patchwork::wrap_plots(pl, ncol = facet_ncol) +
        patchwork::plot_annotation(
            title = title,
            subtitle = subtitle
        ) & theme_default_survex()
}
