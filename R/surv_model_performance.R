#' Helper functions for `model_performance.R`
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param ... other parameters, currently ignored
#' @param times a numeric vector, time points at which ROC curves are calculated if `type == "roc"` or at which metrics are calculated if `type == "metrics"`. Note: if `type=="roc"` this parameter is obligatory
#' @param type character, either `"metrics"` which calculates performance metrics or `"roc"` which calculates ROC curves at specific time points
#' @param metrics a named vector containing the metrics to be calculated. The values should be standardized loss functions. The functions can be supplied manually but has to have these named parameters (`y_true`, `risk`, `surv`, `times`), where `y_true` represents the `survival::Surv` object with observed times and statuses, `risk` is the risk score calculated by the model, and `surv` is the survival function for each observation evaluated at `times`.
#'
#' @return Either a list when `type == "metrics"` or a data.frame if `type == "roc"`
#'
#' @keywords internal
surv_model_performance <- function(explainer, ..., times = NULL, type = "metrics", metrics = NULL) {
    newdata <- explainer$data
    if (type == "metrics") {
    if (is.null(times)) times <- explainer$times
    sf <- explainer$predict_survival_function(explainer$model, newdata, times)
    risk <- explainer$predict_function(explainer$model, newdata)
    y <- explainer$y
    ret_list <- lapply(metrics, function(x) { output <- x(y, risk, sf, times)
                                                  attr(output, "loss_type") <- attr(x, "loss_type")
                                                  output})

    ret_list <- list(result = ret_list, eval_times = times)

    class(ret_list) <- c("surv_model_performance", class(ret_list))
    attr(ret_list, "label") <- explainer$label

    ret_list
    }

    else {
        if (is.null(times)) stop("Times cannot be NULL for type `roc`")
        rocs <- lapply(times, function(time) {
            censored_earlier_mask <- (explainer$y[, 1] < time & explainer$y[, 2] == 0)
            event_later_mask <- explainer$y[, 1] > time
            newdata_t <- newdata[!censored_earlier_mask, ]
            labels <- explainer$y[,2]
            labels[event_later_mask] <- 0
            labels <- labels[!censored_earlier_mask]
            scores <- explainer$predict_survival_function(explainer$model, newdata_t, time)
            labels <- labels[order(scores, decreasing = FALSE)]
            cbind(time = time, data.frame(TPR = cumsum(labels) / sum(labels),
                                          FPR = cumsum(!labels) / sum(!labels), labels))
        })


        rocs_df <- do.call(rbind, rocs)

        class(rocs_df) <- c("surv_model_performance_rocs", class(rocs_df))
        attr(rocs_df, "label") <- explainer$label
        rocs_df
    }


}
