#' Plot 2-Dimensional Model Profile for Survival Models
#'
#' This function plots objects of class `"model_profile_2d_survival"` created
#' using the `model_profile_2d()` function.
#'
#' @param x an object of class `model_profile_2d_survival` to be plotted
#' @param ... additional objects of class `model_profile_2d_survival` to be plotted together
#' @param variables list of character vectors of length 2, names of pairs of variables to be plotted
#' @param times numeric vector, times for which the profile should be plotted, the times must be present in the 'times' field of the explainer. If `NULL` (default) then the median time from the explainer object is used.
#' @param facet_ncol number of columns for arranging subplots
#' @param title character, title of the plot. `'default'` automatically generates either "2D partial dependence survival profiles" or "2D accumulated local effects survival profiles" depending on the explanation type.
#' @param subtitle  character, subtitle of the plot, `'default'` automatically generates "created for the XXX model", where XXX is the explainer labels, if `marginalize_over_time = FALSE`, time is also added to the subtitle
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#'
#' @return A collection of `ggplot` objects arranged with the `patchwork` package.
#'
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' cph_model_profile_2d <- model_profile_2d(cph_exp,
#'                                         variables = list(c("age", "celltype"),
#'                                                          c("age", "karno")))
#' head(cph_model_profile_2d$result)
#' plot(cph_model_profile_2d, variables = list(c("age", "celltype")), times = 88.5)
#'
#' cph_model_profile_2d_ale <- model_profile_2d(cph_exp,
#'                                         variables = list(c("age", "karno")),
#'                                         type = "accumulated")
#' head(cph_model_profile_2d_ale$result)
#' plot(cph_model_profile_2d_ale, times = c(4, 88.5), marginalize_over_time = TRUE)
#' }
#'
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

    if (!is.null(variables)) {
        variables <- intersect(x$variables, variables)
        if (length(variables) == 0)
            stop(paste0(
                "variables do not overlap with ",
                paste(x$variables, collapse = ", ")
            ))
    } else {
        variables <- x$variables
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
        labels[[i]] <- c(this_title, rep("", length(variables)-1))
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
        warning("Plot will be prepared for the median time point from the `times` vector. For another time point, set the value of `times`.")
    }

    if (!marginalize_over_time && length(times) > 1) {
        times <- times[1]
        warning("Plot will be prepared for the first time point in the `times` vector. For aggregation over time, set the option `marginalize_over_time = TRUE`.")
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
    df_time$`_times_` <- NULL
    if (marginalize_over_time){
        df_time <- aggregate(`_yhat_`~., data=df_time, FUN=mean)
    }
    sf_range <- range(df_time$`_yhat_`)

    pl <- lapply(seq_along(variables), function(i){
        variable_pair <- variables[[i]]
        df <- df_time[df_time$`_v1name_` == variable_pair[1] &
                          df_time$`_v2name_` == variable_pair[2],]
        if (any(df$`_v1type_` == "numerical"))
            df$`_v1value_` <- as.numeric(as.character(df$`_v1value_`))
        else if (any(df$`_v1type_` == "categorical"))
            df$`_v1value_` <- as.character(df$`_v1value_`)
        if (any(df$`_v2type_` == "numerical"))
            df$`_v2value_` <- as.numeric(as.character(df$`_v2value_`))
        else if (any(df$`_v2type_` == "categorical"))
            df$`_v2value_` <- as.character(df$`_v2value_`)
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




