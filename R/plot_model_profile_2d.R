#' @export
plot.model_profile_2d_survival <- function(x,
                                           ...,
                                           variables = NULL,
                                           times = NULL,
                                           marginalize_over_time = FALSE,
                                           facet_ncol = NULL,
                                           title = "default",
                                           subtitle = "default",
                                           colors = NULL){

    explanations_list <- c(list(x), list(...))
    num_models <- length(explanations_list)
    if (title == "default"){
        if (x$type == "partial")
            title <- "2D partial dependence survival profiles"
        if (x$type == "accumulated")
            title <- "2D accumulated local effects survival profiles"
    }

    all_variables <- x$variables
    if (!is.null(variables)) {
        all_variables <- intersect(all_variables, variables)
        if (length(all_variables) == 0)
            stop(paste0(
                "variables do not overlap with ",
                paste(all_variables, collapse = ", ")
            ))
    }

    if (is.null(colors))
        colors <- c("#c7f5bf", "#8bdcbe", "#46bac2", "#4378bf", "#371ea3")

    if (num_models == 1){
        result <- prepare_model_profile_2d_plots(x,
                                                variables = variables,
                                                times = times,
                                                marginalize_over_time = marginalize_over_time,
                                                facet_ncol = facet_ncol,
                                                title = title,
                                                subtitle = subtitle,
                                                colors = colors
                                                )
        return(result)
    }

    return_list <- list()
    labels <- list()
    for (i in 1:num_models){
        this_title <- unique(explanations_list[[i]]$result$`_label_`)
        return_list[[i]] <-  prepare_model_profile_2d_plots(explanations_list[[i]],
                                                          variables = variables,
                                                          times = times,
                                                          marginalize_over_time = marginalize_over_time,
                                                          facet_ncol = 1,
                                                          title = this_title,
                                                          subtitle = subtitle,
                                                          colors = colors)
        labels[[i]] <- c(this_title, rep("", length(return_list[[i]]$patches)-2))
    }

    labels <- unlist(labels)
    patchwork::wrap_plots(return_list, nrow = 1, tag_level="keep") +
        patchwork::plot_annotation(title, tag_levels = list(labels)) & theme_default_survex()
}


prepare_model_profile_2d_plots <- function(x,
                                           variables,
                                           times,
                                           marginalize_over_time,
                                           facet_ncol,
                                           title,
                                           subtitle,
                                           colors
){
    if (is.null(times)) {
        times <- quantile(x$eval_times, p = 0.5, type = 1)
    }

    if (!marginalize_over_time) {
        times <- times[1]
    }

    if (!all(times %in% x$eval_times)) {
        stop(paste0(
            "For one of the provided times the explanations has not been calculated or found.
         Please modify the times argument in your explainer or use only values from the following:  ",
            paste(x$eval_times, collapse = ", ")
        ))
    }

    all_profiles <- x$result
    df_time <- all_profiles[all_profiles$`_times_` %in% times, ]

    sf_range <- range(df_time$`_yhat_`)

    pl <- lapply(seq_along(x$variables), function(i){
        variable_pair <- x$variables[[i]]
        df <- df_time[df_time$`_v1name_` == variable_pair[1] &
                          df_time$`_v2name_` == variable_pair[2],]
        if (any(df$`_v1type_` == "numerical"))
            df$`_v1value_` <- as.numeric(as.character(df$`_v1value_`))
        if (any(df$`_v2type_` == "numerical"))
            df$`_v2value_` <- as.numeric(as.character(df$`_v2value_`))
        xlabel <- unique(df$`_v1name_`)
        ylabel <- unique(df$`_v2name_`)

        if (x$type == "partial"){
            p <- with(df, {
                ggplot(df,
                       aes(x = `_v1value_`, y = `_v2value_`, fill = `_yhat_`)) +
                    geom_tile() +
                    scale_fill_gradientn(name = "PDP value",
                                         colors = rev(grDevices::colorRampPalette(colors)(10)),
                                         limits = sf_range) +
                    labs(x = xlabel, y = ylabel) +
                    theme(legend.position = "top") +
                    facet_wrap(~paste(`_v1name_`, `_v2name_`, sep = " : "))
            })
        } else {
            p <-  with(df, {
                ggplot(df, aes(x = `_v1value_`, y = `_v2value_`, fill = `_yhat_`)) +
                    geom_rect(aes(ymin = `_bottom_`, ymax = `_top_`,
                                  xmin = `_left_`, xmax = `_right_`)) +
                    scale_fill_gradientn(name = "ALE value",
                                         colors = rev(grDevices::colorRampPalette(colors)(10)),
                                         limits = sf_range) +
                    labs(x = xlabel, y = ylabel) +
                    theme(legend.position = "top") +
                    facet_wrap(~paste(`_v1name_`, `_v2name_`, sep = " : "))
                })
        }

        if (i != length(x$variables))
            p <- p + guides(fill = "none")
        return(p)
    })
    if (!is.null(subtitle) && subtitle == "default") {
        labels <-
            paste0(unique(all_profiles$`_label_`), collapse = ", ")
        subtitle <- paste0("created for the ", labels, " model")
        if (!marginalize_over_time)
            subtitle <- paste0(subtitle, " and t = ", times)
    }

    patchwork::wrap_plots(pl, ncol = facet_ncol) &
        patchwork::plot_annotation(title = title,
                                   subtitle = subtitle) & theme_default_survex() &
        plot_layout(guides = "collect")
}




