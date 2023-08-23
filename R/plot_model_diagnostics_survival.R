#' Plot Model Diagnostics for Survival Models
#'
#' @export
plot.model_diagnostics_survival <- function(x,
                                            ...,
                                            type = "deviance",
                                            colors = "red"){
    lapply(list(x, ...), function(x) {
        if (!inherits(x, "model_diagnostics_survival")) {
            stop("All ... must be objects of class `model_diagnostics_survival`.")
        }
    })
    explanations_list <- c(list(x), list(...))
    result_list <- lapply(explanations_list, function(x) x$result)
    df <- do.call(rbind, result_list)
    print(df)
    ggplot(df, aes(x = time, y = deviance_residuals, color = status)) +
        geom_point() +
        theme_default_survex() +
        facet_wrap(~label)
}
