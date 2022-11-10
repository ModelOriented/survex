#' Plot SurvLIME Explanations for Survival Models
#'
#' This functions plots objects of class `surv_lime` - LIME explanations of survival models
#' created using `predict_parts(..., type="survlime")` function.
#'
#' @param x an object of class `"surv_lime"` to be plotted
#' @param type character, either "coefficients" or "local_importance" (default), selects the type of plot
#' @param show_survival_function logical, if the survival function of the explanations should be plotted next to the barplot
#' @param ... other parameters currently ignored
#' @param title character, title of the plot
#' @param subtitle character, subtitle of the plot, `'default'` automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
#'
#' @return An object of the class `ggplot`.
#'
#' @family functions for plotting 'predict_parts_survival' objects
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#' model <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' exp <- explain(model)
#'
#' p_parts_lime <- predict_parts(exp, veteran[1, -c(3, 4)], type = "survlime")
#' plot(p_parts_lime)
#'
#' @export
plot.surv_lime <- function(x,
                           type = "local_importance",
                           show_survival_function = TRUE,
                           ...,
                           title = "SurvLIME",
                           subtitle = "default",
                           colors = NULL) {
    if (!type %in% c("coefficients", "local_importance"))
        stop("Type should be one of `coefficients`, `local_importance`")

    x$beta <- x$result


    df <- data.frame(variable_names = names(x$variable_values),
                     variable_values = x$variable_values,
                     beta = x$beta,
                     sign_beta = as.factor(sign(x$result)),
                     sign_local_importance = as.factor(sign(x$beta * x$variable_values)),
                     local_importance = x$beta * x$variable_values)

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", attr(x, "label"), " model")
    }

    if (type == "coefficients") {
        x_lab <- "SurvLIME coefficients"
        y_lab <- ""
        pl <- with(df, {
        ggplot(data = df, aes(x = beta, y = reorder(variable_names, beta, abs), fill = sign_beta)) +
            geom_col() +
            scale_fill_manual("", values = c("#f05a71", "#8bdcbe"))
        })
    }


    if (type == "local_importance") {
        x_lab <- "SurvLIME local importance"
        y_lab <- ""
        pl <- ggplot(data = df, aes_string(x = "local_importance", y = "reorder(variable_names, local_importance, abs)", fill = "sign_local_importance")) +
            geom_col() +
            scale_fill_manual("", values = c("#f05a71", "#ffffff", "#8bdcbe"))
    }
    pl <- pl + theme_drwhy_vertical() +
        labs(title = title, subtitle = subtitle) +
        xlab(x_lab) +
        ylab(y_lab) +
        theme(legend.position = "none")

    sf_df <- data.frame(times = c(x$black_box_sf_times, x$expl_sf_times),
                        sfs = c(x$black_box_sf, x$expl_sf),
                        type = c(rep("black box survival function", length(x$black_box_sf)), rep("SurvLIME explanation survival function", length(x$expl_sf))))
    if (show_survival_function) {
        pl2 <- ggplot(data = sf_df, aes_string(x = "times", y = "sfs", group = "type", color = "type")) +
            geom_line(size = 0.8) +
            theme_drwhy() +
            xlab("") +
            ylab("survival function value") +
            scale_color_manual("", values = generate_discrete_color_scale(2, colors))
        return(patchwork::wrap_plots(pl, pl2, nrow = 1, widths = c(3, 5)))
    }


}