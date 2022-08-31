#' Helper functions for `model_performance.R`
#'
#' @param explainer a model to be explained, preprocessed by the `explain` function
#' @param ... other parameters, currently ignored
#' @param times a numeric vector, time points at which ROC curves are calculated if `type == "roc"` or at which metrics are calculated if `type == "metrics"`. Note: if `type=="roc"` this parameter is obligatory
#' @param type character, either `"metrics"` which calculates performance metrics or `"roc"` which calculates ROC curves at specific time points
#'
#' @return Either a list when `type == "metrics"` or a data.frame if `type == "roc"`
#'
#' @keywords internal
surv_model_performance <- function(explainer, ..., times = NULL, type = "metrics") {
    newdata <- explainer$data
    if (type == "metrics") {
    if (is.null(times)) times <- explainer$times
    sf <- explainer$predict_survival_function(explainer$model, newdata, times)
    risk <- explainer$predict_function(explainer$model, newdata)
    y <- explainer$y


    brier_score <- brier_score(y, risk, sf, times)
    auc <- cd_auc(y, risk, sf,  times)
    cindex <- c_index(y, risk)
    iauc <- integrated_cd_auc(auc = auc, times = times)
    ibs <- integrated_brier_score(times = times, brier = brier_score)

    ret_list <- list(
        eval_times = times,
        brier_score = brier_score,
        auc = auc,
        cindex = cindex,
        iauc = iauc,
        integrated_brier_score = ibs
    )

    class(ret_list) <- c("surv_model_performance", class(ret_list))
    attr(ret_list, "label") <- explainer$label

    ret_list
    }

    else {
        if (is.null(times)) stop("Times cannot be NULL for type `roc`")
        rocs <- lapply(times, function(time) {
            labels <- 1 - explainer$y[, 2]
            scores <- explainer$predict_survival_function(explainer$model, newdata, time)
            labels <- labels[order(scores, decreasing = TRUE)]
            cbind(time = time, data.frame(TPR = cumsum(labels) / sum(labels), FPR = cumsum(!labels) / sum(!labels), labels))
        })

        rocs_df <- do.call(rbind, rocs)

        class(rocs_df) <- c("surv_model_performance_rocs", class(rocs_df))
        attr(rocs_df, "label") <- explainer$label
        rocs_df
    }


}
