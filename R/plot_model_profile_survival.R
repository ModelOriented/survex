#' Plot Model Profile for Survival Models
#'
#' This function plots objects of class `"model_profile_survival"` created
#' using the `model_profile()` function.
#'
#' @param x an object of class `model_profile_survival` to be plotted
#' @param ... additional objects of class `"model_profile_survival"` to be plotted together
#' @param variables character, names of the variables to be plotted
#' @param variable_type character, either `"numerical"`, `"categorical"` or `NULL` (default), select only one type of variable for plotting, or leave `NULL` for all
#' @param facet_ncol number of columns for arranging subplots
#' @param numerical_plot_type character, either `"lines"`, or `"contours"` selects the type of numerical variable plots
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' @param rug character, one of `"all"`, `"events"`, `"censors"`, `"none"` or `NULL`. Which times to mark on the x axis in `geom_rug()`.
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
#' }
#'
#' @export
plot.model_profile_survival <- function(x,
                                        ...,
                                        variables = NULL,
                                        variable_type = NULL,
                                        facet_ncol = NULL,
                                        numerical_plot_type = "lines",
                                        title = "Partial dependence survival profile",
                                        subtitle = "default",
                                        colors = NULL,
                                        rug = "all",
                                        rug_colors = c("#dd0000", "#222222")) {
    explanations_list <- c(list(x), list(...))
    num_models <- length(explanations_list)

    if (num_models == 1) {
        result <- prepare_model_profile_plots(x,
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
        return_list[[i]] <- prepare_model_profile_plots(explanations_list[[i]],
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


#' @export
plot2 <- function(x, ...) UseMethod("plot2")

#' @export
plot2.model_profile_survival <- function(x,
                                         variable,
                                         times = NULL,
                                         marginalize_over_time = FALSE,
                                         plot_type = "pdp+ice",
                                         ...,
                                         title = "Partial dependence profile",
                                         subtitle = "default",
                                         colors = NULL) {
    if (is.null(variable) || !is.character(variable)) {
        stop("A variable must be specified by name")
    }

    if (length(variable) > 1) {
        stop("Only one variable can be specified")
    }

    if (!variable %in% x$result$`_vname_`) {
        stop(paste0("Variable ", variable, " not found"))
    }

    if (is.null(times)) {
        times <- quantile(x$eval_times, p = 0.5, type = 1)
    }

    if (!all(times %in% x$eval_times)) {
        stop(paste0(
            "For one of the provided times the explanations has not been calculated not found.
         Please modify the times argument in your explainer or use only values from the following:  ",
            paste(x$eval_times, collapse = ", ")
        ))
    }

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", unique(x$result$`_label_`), " model")
    }

    # Select relevant information from the pdp result
    pdp_df <- x$result[(x$result$`_vname_` == variable) & (x$result$`_times_` %in% times), c("_x_", "_yhat_")]
    colnames(pdp_df) <- c(variable, "pd")

    # Select relevant information from the ceteris paribus profiles
    # TODO: REMOVE THIS WHEN ID IS FIXED
    # for (i in 1:length(x$cp_profiles)) {
    #     return(unique(x$cp_profiles$result$`_vname_`))
    # }

    ice_df <- x$cp_profiles$result[(x$cp_profiles$result$`_vname_` == variable) &
        (x$cp_profiles$result$`_times_` %in% times), ]
    ice_df$`_times_` <- NULL
    ice_df$`_vname_` <- NULL
    ice_df$`_vtype_` <- NULL
    ice_df$`_label_` <- NULL

    colnames(ice_df)[colnames(ice_df) == "_ids_"] <- "id"
    colnames(ice_df)[colnames(ice_df) == "_yhat_"] <- "predictions"


    feature_name_sym <- sym(variable)
    y_floor_pd <- floor(min(pdp_df[, "pd"]) * 10) / 10
    y_ceiling_pd <- ceiling(max(pdp_df[, "pd"]) * 10) / 10

    single_timepoint <- ((length(times) == 1) || marginalize_over_time)

    return(ice_df)

    if (unique(x$result[x$result$`_vname_` == variable, "_vtype_"]) == "categorical") {
        pl <- plot_pdp_cat(
            pdp_dt = pdp_df,
            ice_dt = ice_df,
            data_dt = NULL,
            feature_name_sym,
            y_floor_ice = NULL,
            y_ceiling_ice = NULL,
            y_floor_pd = y_floor_pd,
            y_ceiling_pd = y_ceiling_pd,
            plot_type = plot_type,
            single_timepoint = single_timepoint
        )
    } else {
        pdp_df[, 1] <- as.numeric(as.character(pdp_df[, 1]))
        pl <- plot_pdp_num(
            pdp_dt = pdp_df,
            ice_dt = ice_df,
            data_dt = NULL,
            feature_name_sym,
            y_floor_ice = NULL,
            y_ceiling_ice = NULL,
            y_floor_pd = y_floor_pd,
            y_ceiling_pd = y_ceiling_pd,
            plot_type = plot_type,
            single_timepoint = single_timepoint
        )
    }
    pl
}

plot_pdp_num <- function(pdp_dt,
                         ice_dt,
                         data_dt,
                         feature_name_sym,
                         y_floor_ice,
                         y_ceiling_ice,
                         y_floor_pd,
                         y_ceiling_pd,
                         plot_type,
                         single_timepoint) {
    if (single_timepoint == TRUE) { ## single timepoint
        if (plot_type == "ice") {
            ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                geom_line(alpha = 0.2, mapping = aes(group = id)) +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = "jitter")
        }
        # PDP + ICE
        else if (plot_type == "pdp+ice") {
            ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                geom_line(mapping = aes(group = id), alpha = 0.2) +
                geom_line(data = pdp_dt, aes(x = !!feature_name_sym, y = pd), linewidth = 2, color = "gold") #+
            # geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = "jitter")
        }
        # PDP
        else if (plot_type == "pdp") {
            ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
                geom_line() +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_pd, y_ceiling_pd)
        }
    } else { ## multiple timepoints
        if (plot_type == "ice") {
            ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                geom_line(alpha = 0.2, mapping = aes(group = interaction(id, time), color = time)) +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides = "b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_ice, y_ceiling_ice)
        }
        # PDP + ICE
        else if (plot_type == "pdp+ice") {
            ggplot() +
                geom_line(data = ice_dt, aes(x = !!feature_name_sym, y = predictions, group = interaction(id, time), color = time), alpha = 0.1) +
                geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, color = time), linewidth = 1.5, lineend = "round", linejoin = "round") +
                geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, group = time), color = "black", linewidth = 0.5, linetype = "dashed", lineend = "round", linejoin = "round") +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides = "b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_ice, y_ceiling_ice)
        }
        # PDP
        else if (plot_type == "pdp") {
            ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
                geom_line(aes(color = time)) +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_pd, y_ceiling_pd)
        }
    }
}

plot_pdp_cat <- function() {

}




prepare_model_profile_plots <- function(x,
                                        variables = NULL,
                                        variable_type = NULL,
                                        facet_ncol = NULL,
                                        numerical_plot_type = "lines",
                                        title = "Partial dependence survival profile",
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
        labels <-
            paste0(unique(aggregated_profiles$`_label_`), collapse = ", ")
        subtitle <- paste0("created for the ", labels, " model")
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

    pl <- plot_individual_ceteris_paribus_survival(aggregated_profiles, variables, colors, numerical_plot_type, rug_df, rug, rug_colors)

    patchwork::wrap_plots(pl, ncol = facet_ncol) +
        patchwork::plot_annotation(
            title = title,
            subtitle = subtitle
        ) & theme_default_survex()
}
