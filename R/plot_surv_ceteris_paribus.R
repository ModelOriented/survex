prepare_ceteris_paribus_plots <- function(x,
                                          colors = NULL,
                                          variable_type = NULL,
                                          facet_ncol = NULL,
                                          variables = NULL,
                                          numerical_plot_type = "lines",
                                          title = "default",
                                          subtitle = "default",
                                          rug = "all",
                                          rug_colors = c("#dd0000", "#222222")) {
    rug_df <- data.frame(times = x$event_times, statuses = as.character(x$event_statuses), label = unique(x$result$`_label_`))
    obs <- as.data.frame(x$variable_values)
    center <- x$center
    x <- x$result

    all_profiles <- x
    class(all_profiles) <- "data.frame"

    all_profiles$`_ids_` <- factor(all_profiles$`_ids_`)

    # extract labels to use in the default subtitle
    if (!is.null(subtitle) && subtitle == "default") {
        labels <- unique(all_profiles$`_label_`)
        endword <- ifelse(length(labels) > 1, " models", " model")
        subtitle <- paste0("created for the ", paste0(labels, collapse = ", "), endword)
    }

    # variables to use
    all_variables <- na.omit(as.character(unique(all_profiles$`_vname_`)))
    if (!is.null(variables)) {
        variables <- intersect(all_variables, variables)
        if (length(variables) == 0) {
            stop(paste0(
                "variables do not overlap with ",
                paste(all_variables, collapse = ", ")
            ))
        }
        all_variables <- variables
    }


    if (!is.null(variable_type) && variable_type == "numerical") {
        x <- x[x$`_vtype_` == "numerical", ]
    }

    if (!is.null(variable_type) && variable_type == "categorical") {
        x <- x[x$`_vtype_` == "categorical", ]
    }

    all_variables <- intersect(all_variables, unique(x$`_vname_`))

    lsc <- lapply(all_variables, function(sv) {
        tmp <- x[
            x$`_vname_` == sv,
            c(
                sv,
                "_times_",
                "_vname_",
                "_vtype_",
                "_yhat_",
                "_label_",
                "_ids_"
            )
        ]
        key <- obs[, sv, drop = FALSE]

        tmp$`_real_point_` <- tmp[, sv] == key[, sv]

        colnames(tmp)[1] <- "_x_"
        tmp$`_x_` <- as.character(tmp$`_x_`)
        tmp
    })

    ceteris_paribus_with_observation <- do.call(rbind, lsc)

    pl <- plot_individual_ceteris_paribus_survival(
        all_profiles = ceteris_paribus_with_observation,
        variables = all_variables,
        colors = colors,
        numerical_plot_type = numerical_plot_type,
        rug_df = rug_df,
        rug = rug,
        rug_colors = rug_colors,
        center = center
    )

    patchwork::wrap_plots(pl, ncol = facet_ncol) +
        patchwork::plot_annotation(
            title = title,
            subtitle = subtitle
        ) & theme_default_survex()
}

#' @import ggplot2
#' @importFrom graphics par
plot_individual_ceteris_paribus_survival <- function(all_profiles,
                                                     variables,
                                                     colors,
                                                     numerical_plot_type,
                                                     rug_df,
                                                     rug,
                                                     rug_colors,
                                                     center) {
    pl <- lapply(variables, function(var) {
        df <- all_profiles[all_profiles$`_vname_` == var, ]

        if (unique(df$`_vtype_`) == "numerical") {
            if (!is.null(colors)) {
                scale_cont <- colors
            } else {
                scale_cont <- c(
                    low = "#9fe5bd",
                    mid = "#46bac2",
                    high = "#371ea3"
                )
            }
            if (numerical_plot_type == "lines") {
                base_plot <- with(df, {
                    ggplot(
                        df,
                        aes(
                            x = `_times_`,
                            y = `_yhat_`,
                            group = `_x_`,
                            color = as.numeric(as.character(`_x_`))
                        )
                    ) +
                        geom_line() +
                        scale_colour_gradient2(
                            name = paste0(unique(df$`_vname_`), " value"),
                            low = scale_cont[1],
                            high = scale_cont[3],
                            mid = scale_cont[2],
                            midpoint = median(as.numeric(as.character(df$`_x_`)))
                        ) +
                        geom_line(
                            data = df[df$`_real_point_`, ], color =
                                "red", linewidth = 0.8
                        ) +
                        labs(x = "time", y = "centered profile value") +
                        xlim(c(0, NA)) +
                        theme_default_survex() +
                        facet_wrap(~`_vname_`)
                })
                if (!center) {
                    base_plot <- base_plot + ylim(c(0, 1)) + ylab("survival function value")
                }
            } else {
                base_plot <- with(df, {
                    ggplot(
                        df,
                        aes(
                            x = `_times_`,
                            y = as.numeric(as.character(`_x_`)),
                            z = `_yhat_`
                        )
                    )
                })
                if (!center) {
                    base_plot <- base_plot +
                        geom_contour_filled(binwidth = 0.1) +
                        scale_fill_manual(
                            name = "SF value",
                            values = grDevices::colorRampPalette(c("#c7f5bf", "#8bdcbe", "#46bac2", "#4378bf", "#371ea3"))(10),
                            drop = FALSE
                        ) +
                        guides(fill = guide_colorsteps(
                            direction = "horizontal",
                            barwidth = 0.5 * unit(par("pin")[1], "in"),
                            barheight = 0.02 * unit(par("pin")[2], "in"),
                            reverse = TRUE,
                            show.limits = TRUE
                        )) +
                        labs(x = "time", y = "variable value") +
                        xlim(c(0, NA)) +
                        theme_default_survex() +
                        facet_wrap(~`_vname_`)
                } else {
                    base_plot <- base_plot +
                        geom_contour_filled(bins = 10) +
                        scale_fill_manual(
                            name = "centered profile value",
                            values = grDevices::colorRampPalette(c("#c7f5bf", "#8bdcbe", "#46bac2", "#4378bf", "#371ea3"))(10),
                            drop = FALSE
                        ) +
                        guides(fill = guide_colorsteps(
                            direction = "horizontal",
                            barwidth = 0.5 * unit(par("pin")[1], "in"),
                            barheight = 0.02 * unit(par("pin")[2], "in"),
                            reverse = TRUE,
                            show.limits = TRUE,
                            label.theme = element_text(size = 7)
                        )) +
                        labs(x = "time", y = "variable value") +
                        xlim(c(0, NA)) +
                        theme_default_survex() +
                        facet_wrap(~`_vname_`)
                }
                if (any(df$`_real_point_`)) {
                    range_time <- range(df["_times_"])
                    var_val <- as.numeric(unique(df[df$`_real_point_`, "_x_"]))
                    base_plot <- base_plot + geom_segment(aes(x = range_time[1], y = var_val, xend = range_time[2], yend = var_val), color = "red")
                }
                base_plot
            }
        } else {
            n_colors <- length(unique(df$`_x_`))
            base_plot <- with(df, {
                ggplot(
                    df,
                    aes(
                        x = `_times_`,
                        y = `_yhat_`,
                        group = `_x_`,
                        color = `_x_`
                    )
                ) +
                    geom_line(
                        data = df[!df$`_real_point_`, ],
                        linewidth = 0.8
                    ) +
                    geom_line(
                        data = df[df$`_real_point_`, ],
                        linewidth = 0.8, linetype = "longdash"
                    ) +
                    scale_color_manual(
                        name = paste0(unique(df$`_vname_`), " value"),
                        values = generate_discrete_color_scale(n_colors, colors)
                    ) +
                    theme_default_survex() +
                    labs(x = "time", y = "centered profile value") +
                    xlim(c(0, NA)) +
                    facet_wrap(~`_vname_`)
            })
            if (!center) {
                base_plot <- base_plot + ylim(c(0, 1)) + ylab("survival function value")
            }
        }

        return_plot <- add_rug_to_plot(base_plot, rug_df, rug, rug_colors)

        return(return_plot)
    })

    pl
}
