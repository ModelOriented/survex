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
#' @param max_vars maximum number of variables to be plotted (least important variables are ignored)
#' @param colors character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
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
#' p_parts_lime <- predict_parts(exp, veteran[1, -c(3, 4)], type = "survlime")
#' plot(p_parts_lime)
#' }
#' @export
plot.surv_lime <- function(x,
                           type = "local_importance",
                           show_survival_function = TRUE,
                           ...,
                           title = "SurvLIME",
                           subtitle = "default",
                           max_vars = 7,
                           colors = NULL) {
    if (!type %in% c("coefficients", "local_importance"))
        stop("Type should be one of `coefficients`, `local_importance`")

    local_importance <- as.numeric(x$result) * as.numeric(x$variable_values)
    df <- data.frame(variable_names = names(x$variable_values),
                     variable_values = as.numeric(x$variable_values),
                     beta = as.numeric(x$result),
                     sign_beta = as.factor(sign(as.numeric(x$result))),
                     sign_local_importance = as.factor(sign(local_importance)),
                     local_importance = local_importance)

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", attr(x, "label"), " model")
    }

    if (type == "coefficients") {
        x_lab <- "SurvLIME coefficients"
        y_lab <- ""
        df <- df[head(order(abs(df$beta), decreasing=TRUE), max_vars),]
        pl <- with(df, {
        ggplot(data = df, aes(x = beta, y = reorder(variable_names, beta, abs), fill = sign_beta)) +
            geom_col() +
            scale_fill_manual("", values = c("-1"="#f05a71", "0"="#ffffff", "1"="#8bdcbe"))
        })
    }


    if (type == "local_importance") {
        x_lab <- "SurvLIME local importance"
        y_lab <- ""
        df <- df[head(order(abs(df$local_importance), decreasing=TRUE), max_vars),]
        pl <- with(df,{

            ggplot(data = df, aes(x = local_importance, y = reorder(variable_names, local_importance, abs), fill = sign_local_importance)) +
            geom_col() +
            scale_fill_manual("", values = c("-1"="#f05a71", "0"="#ffffff", "1"="#8bdcbe"))

        })
    }
    pl <- pl + theme_vertical_default_survex() +
        labs(title = title, subtitle = subtitle) +
        xlab(x_lab) +
        ylab(y_lab) +
        theme(legend.position = "none")

    sf_df <- data.frame(times = c(x$black_box_sf_times, x$expl_sf_times),
                        sfs = c(x$black_box_sf, x$expl_sf),
                        type = c(rep("black box survival function", length(x$black_box_sf)), rep("SurvLIME explanation survival function", length(x$expl_sf))))
    if (show_survival_function) {
        pl2 <- with(sf_df,{

            ggplot(data = sf_df, aes(x = times, y = sfs, group = type, color = type)) +
            geom_line(linewidth = 0.8, size = 0.8) +
            theme_default_survex() +
            xlab("") +
            xlim(c(0,NA))+
            ylab("survival function value") +
            scale_color_manual("", values = generate_discrete_color_scale(2, colors))
        })
        return(patchwork::wrap_plots(pl, pl2, nrow = 1, widths = c(3, 5)))
    }


}
