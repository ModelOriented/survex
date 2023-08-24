#' Plot Model Diagnostics for Survival Models
#'
#' @export
plot.model_diagnostics_survival <- function(x,
                                            ...,
                                            type = "deviance",
                                            xvariable = "index",
                                            smooth = as.logical(xvariable != "index"),
                                            title = "Model diagnostics",
                                            subtitle = "default",
                                            facet_ncol =NULL,
                                            colors = c("#160e3b", "#f05a71", "#ceced9")){
    lapply(list(x, ...), function(x) {
        if (!inherits(x, "model_diagnostics_survival")) {
            stop("All ... must be objects of class `model_diagnostics_survival`.")
        }
    })
    explanations_list <- c(list(x), list(...))
    result_list <- lapply(explanations_list, function(x) x$result)
    n_observations <- print(sapply(result_list, nrow))
    df <- do.call(rbind, result_list)

    if (!is.null(subtitle) && subtitle == "default") {
        subtitle <- paste0("created for the ", paste(unique(df$label), collapse = ", "),
                          ifelse(length(unique(df$label)) > 1, " models", " model"))
    }

    if (type %in% c("deviance", "martingale")){
        if (!xvariable %in% c("index", colnames(df))){
            stop(paste("`xvariable`", xvariable, "not found"))
        }
        df$y <- switch(type,
                       "deviance" = df$deviance_residuals,
                       "martingale" = df$martingale_residuals)
        df$x <- switch(xvariable,
                       "index" = unlist(lapply(n_observations, function(x) seq_len(x))),
                       df[[xvariable]])

        pl <- ggplot(df, aes(x = x, y = y, color = status)) +
            geom_hline(yintercept = 0, color = colors[3], lty = 2, linewidth = 1) +
            geom_point() +
            theme_default_survex() +
            scale_color_manual(values = colors, labels = c("censored", "event")) +
            facet_wrap(~label, ncol = facet_ncol) +
            labs(title = title,
                 subtitle = subtitle,
                 x = xvariable,
                 y = paste(type, "residuals"))
        if (smooth)
            pl <- pl + geom_smooth(se = FALSE, color = colors[3], alpha = 0.5)
        return(pl)
    } else if (type == "Cox-Snell"){
        split_df <- split(df, df$label)
        df_list <- lapply(split_df, function(df_tmp){
            fit_coxsnell <- survival::survfit(survival::Surv(cox_snell_residuals, as.numeric(status)) ~ 1, data=df_tmp)
            confint <- survival:::survfit_confint(fit_coxsnell$cumhaz, fit_coxsnell$std.chaz,
                                       logse=FALSE, "log", 0.95, ulimit = FALSE)
            data.frame(
                "time" = fit_coxsnell$time,
                "cumhaz" = fit_coxsnell$cumhaz,
                "lower" = confint$lower,
                "upper" = confint$upper)
        })

        df <- do.call(rbind, df_list)
        df$label <- sapply(strsplit(rownames(df), "[.]"), function(x) x[1])
        rownames(df) <- NULL

        ggplot(df, aes(x = time, y = cumhaz)) +
            geom_step(color = colors[1], linewidth = 1) +
            geom_step(aes(y = lower), linetype = "dashed", color = colors[1], alpha = 0.8) +
            geom_step(aes(y = upper), linetype = "dashed", color = colors[1], alpha = 0.8) +
            geom_abline(slope = 1, color = colors[2], linewidth = 1) +
            labs(x = "Cox-Snell residuals (pseudo observed times)",
                 y = "Cumulative hazard at pseudo observed times") +
            theme_default_survex() +
            facet_wrap(~label, ncol = facet_ncol) +
            labs(title = title,
                 subtitle = subtitle)
    } else{
        stop('`type` should be one of `deviance`, `martingale` or `Cox-Snell`')
    }
}
