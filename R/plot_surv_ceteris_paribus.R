#' Plot Predict Profile for Survival Models
#'
#' This function plots objects of class `"predict_profile_survival"` created using
#' the `predict_profile()` function.
#'
#' @param x an object of class `predict_profile_survival` to be plotted
#' @param ... additional parameters, unused, currently ignored
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#' @param variable_type character, either `"numerical"`, `"categorical"` or `NULL` (default), select only one type of variable for plotting, or leave `NULL` for all
#' @param facet_ncol number of columns for arranging subplots
#' @param variables character, names of the variables to be plotted
#' @param numerical_plot_type character, either `"lines"`, or `"contours"` selects the type of numerical variable plots
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, if `NULL` automaticaly generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
#'
#' @return A grid of `ggplot` objects arranged with the `gridExtra::grid.arrange` function.
#'
#' @family functions for plotting 'predict_profile_survival' objects
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
plot.surv_ceteris_paribus <- function(x,
                                      ...,
                                      colors = NULL,
                                      variable_type = NULL,
                                      facet_ncol = NULL,
                                      variables = NULL,
                                      numerical_plot_type = "lines",
                                      title = "Ceteris paribus survival profile",
                                      subtitle = NULL) {
    if (!is.null(variable_type))
        check_variable_type(variable_type)

    check_numerical_plot_type(numerical_plot_type)

    obs <- x$variable_values
    x <- x$result

    all_profiles <- x
    class(all_profiles) <- "data.frame"

    all_profiles$`_ids_` <- factor(all_profiles$`_ids_`)

    # extract labels to use in the default subtitle
    if (is.null(subtitle)) {
        labels <- paste0(unique(all_profiles$`_label_`), collapse = ", ")
        subtitle <- paste0("created for the ", labels, " model")
    }

    # variables to use
    all_variables <-
        na.omit(as.character(unique(all_profiles$`_vname_`)))
    if (!is.null(variables)) {
        all_variables <- intersect(all_variables, variables)
        if (length(all_variables) == 0)
            stop(paste0(
                "variables do not overlap with ",
                paste(all_variables, collapse = ", ")
            ))
    }


    if (!is.null(variable_type) && variable_type == "numerical") {
        x <- x[x$`_vtype_` == "numerical", ]
    }

    if (!is.null(variable_type) && variable_type == "categorical") {
        x <- x[x$`_vtype_` == "categorical", ]
    }

    all_variables <- intersect(all_variables, unique(x$`_vname_`))

    if (is.null(facet_ncol))
        facet_ncol <-
        switch(as.character(length(all_variables)), "1" = 1, "2" = 2, 3)

    lsc <- lapply(all_variables, function(sv) {
        tmp <- x[x$`_vname_` == sv,
                 c(sv,
                   "_times_",
                   "_vname_",
                   "_vtype_",
                   "_yhat_",
                   "_label_",
                   "_ids_")]

        key <- obs[, sv, drop = FALSE]

        tmp$`_real_point_` <-
            tmp[, sv] == key[as.character(tmp$`_ids_`), sv]

        colnames(tmp)[1] <- "_x_"
        tmp$`_x_` <- as.character(tmp$`_x_`)
        tmp
    })

    ceteris_paribus_with_observation <- do.call(rbind, lsc)

    pl <- plot_individual_ceteris_paribus_survival(
        all_profiles = ceteris_paribus_with_observation,
        variables = all_variables,
        facet_ncol = facet_ncol,
        colors = colors,
        numerical_plot_type = numerical_plot_type,
        title = title,
        subtitle = subtitle
    )


    titleGrob <-
        grid::textGrob(
            title,
            just = "left",
            x = 0.05,
            gp = grid::gpar(fontsize = 16, col = "#371ea3")
        )
    subtitleGrob <-
        grid::textGrob(
            subtitle,
            just = "left",
            x = 0.05,
            gp = grid::gpar(fontsize = 11, col = "#371ea3")
        )
    grided <-
        do.call(gridExtra::arrangeGrob, c(pl, list(ncol = facet_ncol)))
    margin <- grid::unit(0.5, "line")
    return(
        gridExtra::grid.arrange(
            titleGrob,
            subtitleGrob,
            grided,
            ncol = 1,
            heights = grid::unit.c(
                grid::grobHeight(titleGrob) + 1.2 * margin,
                grid::grobHeight(subtitleGrob) + margin,
                grid::unit(1, "null")
            )
        )
    )
}

#' @importFrom DALEX theme_drwhy
#' @import ggplot2
plot_individual_ceteris_paribus_survival <- function(all_profiles,
                                                     variables,
                                                     facet_ncol,
                                                     colors,
                                                     numerical_plot_type,
                                                     title,
                                                     subtitle) {
    pl <- lapply(variables, function(var) {
        df <- all_profiles[all_profiles$`_vname_` == var, ]

        if (unique(df$`_vtype_`) == "numerical") {
            if (!is.null(colors))
                scale_cont <- colors
            else
                scale_cont <- c(low = "#9fe5bd",
                                mid = "#46bac2",
                                high = "#371ea3")
            if (numerical_plot_type == "lines") {
                ggplot(
                    df,
                    aes_string(
                        x = "`_times_`",
                        y = "`_yhat_`",
                        group = "`_x_`",
                        color = "as.numeric(as.character(`_x_`))"
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
                    geom_line(data = df[df$`_real_point_`, ], color =
                                  "red", size = 0.8) +
                    xlab("") + ylab("survival function value") + ylim(c(0, 1)) +
                    theme_drwhy() +
                    facet_wrap(~`_vname_`)
            } else {
                plt <- ggplot(
                    df,
                    aes_string(
                        x = "`_times_`",
                        y = "as.numeric(as.character(`_x_`))",
                        z = "`_yhat_`"
                    )
                ) +
                    geom_contour_filled(breaks = seq(1, 0, -0.1)) +
                    scale_fill_manual(name = "SF value", values = grDevices::colorRampPalette(c("#c7f5bf", "#8bdcbe", "#46bac2", "#4378bf", "#371ea3"))(10),
                                    labels = seq(1, 0, -0.1)) +
                    guides(fill = guide_legend(nrow = 1, label.position = "top")) +
                    xlab("") + ylab("variable value") +
                    theme_drwhy() +
                    theme(legend.spacing = grid::unit(0.1, 'line')) +
                    facet_wrap(~`_vname_`)
                if (any(df$`_real_point_`)) {
                    range_time <- range(df["_times_"])
                    var_val <- as.numeric(unique(df[df$`_real_point_`, "_x_"]))
                    plt <- plt + geom_segment(aes(x = range_time[1], y = var_val, xend = range_time[2], yend = var_val), color = "red")
                    }
                plt
                }
        } else {
            n_colors <- length(unique(df$`_x_`))
            pl <-
                ggplot(
                    df,
                    aes_string(
                        x = "`_times_`",
                        y = "`_yhat_`",
                        group = "`_x_`",
                        color = "`_x_`"
                    )
                ) +
                geom_line(data = df[!df$`_real_point_`, ],
                          size = 0.8) +
                geom_line(data = df[df$`_real_point_`, ],
                          size = 0.8, linetype = "longdash") +
                scale_color_manual(name = "",
                                   values = generate_discrete_color_scale(n_colors, colors)) +
                facet_wrap(~`_vname_`, scales = "free_x", ncol = facet_ncol) +
                theme_drwhy() +
                xlab("") + ylab("survival function value") + ylim(c(0, 1))
        }
    })
}


check_variable_type <- function(variable_type) {
    if (!(variable_type %in% c("numerical", "categorical")))
        stop("variable_type needs to be 'numerical' or 'categorical'")
}

check_numerical_plot_type <- function(numerical_plot_type) {
    if (!(numerical_plot_type %in% c("lines", "contours")))
        stop("numerical_plot_type needs to be 'lines' or 'contours'")
}
