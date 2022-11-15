utils::globalVariables(c("PredictionSurv"))
#' Compute the Harrell's Concordance index
#'
#' A function to compute the Harrell's concordance index of a survival model.
#'
#'
#' @param y_true a `survival::Surv` object containing the times and statuses of observations for which the metric will be evaluated
#' @param risk a numeric vector of risk scores corresponding to each observation
#' @param surv ignored, left for compatibility with other metrics
#' @param times ignored, left for compatibility with other metrics
#'
#' @return numeric from 0 to 1, higher values indicate better performance
#'
#' @section References:
#' - \[1\] Harrell, F.E., Jr., et al. ["Regression modelling strategies for improved prognostic prediction."](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4780030207) Statistics in Medicine 3.2 (1984): 143-152.
#'
#' @rdname c_index
#' @seealso [loss_one_minus_c_index()]
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' rotterdam <- survival::rotterdam
#' rotterdam$year <- NULL
#' cox_rotterdam_rec <- coxph(Surv(rtime, recur) ~ .,
#'     data = rotterdam,
#'     model = TRUE, x = TRUE, y = TRUE)
#' coxph_explainer <- explain(cox_rotterdam_rec)
#'
#' risk <- coxph_explainer$predict_function(coxph_explainer$model, coxph_explainer$data)
#' c_index(y_true = coxph_explainer$y, risk = risk)
#' }
#'
#' @export
c_index <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL) {

    n_rows <- length(y_true[, 1])

    yi <- matrix(rep(y_true[, 1], n_rows), ncol = n_rows)
    yj <- t(matrix(rep(y_true[, 1], n_rows), ncol = n_rows))

    dj <- t(matrix(rep(y_true[, 2], n_rows), ncol = n_rows))

    ri <- matrix(rep(risk, n_rows), ncol = n_rows)
    rj <- t(matrix(rep(risk, n_rows), ncol = n_rows))

    top <- sum(ifelse(yj < yi & rj > ri, 1, 0) * dj)
    bot <- sum(ifelse(yj < yi, 1, 0) * dj)

    top / bot

}
attr(c_index, "loss_name") <- "C-index"
attr(c_index, "loss_type") <- "risk-based"


#' Calculate the Concordance index loss
#'
#' This function subtracts the C-index metric from one to obtain a loss function whose lower values indicate better model performance (useful for permutational feature importance)
#'
#' @inheritParams c_index
#'
#' @return numeric from 0 to 1, lower values indicate better performance
#'
#' @section References:
#' - \[1\] Harrell, F.E., Jr., et al. ["Regression modelling strategies for improved prognostic prediction."](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4780030207) Statistics in Medicine 3.2 (1984): 143-152.
#'
#' @rdname loss_one_minus_c_index
#' @seealso [c_index()]
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' rotterdam <- survival::rotterdam
#' rotterdam$year <- NULL
#' cox_rotterdam_rec <- coxph(Surv(rtime, recur) ~ .,
#'     data = rotterdam,
#'     model = TRUE, x = TRUE, y = TRUE)
#' coxph_explainer <- explain(cox_rotterdam_rec)
#'
#' risk <- coxph_explainer$predict_function(coxph_explainer$model, coxph_explainer$data)
#' loss_one_minus_c_index(y_true = coxph_explainer$y, risk = risk)
#' }
#' @export
loss_one_minus_c_index <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL) {
    1 - c_index(y_true = y_true, risk = risk, surv = surv, times = times)
}
attr(loss_one_minus_c_index, "loss_name") <- "One minus C-Index"
attr(loss_one_minus_c_index, "loss_type") <- "risk-based"

#' Calculate Brier score
#'
#' A function for calculating the Brier score for a survival model.
#'
#' Brier score is used to evaluate the performance of a survival model, based on the squared distance between the predicted survival function and the actual event time, weighted to account for censored observations.
#'
#'
#' @param y_true a `survival::Surv` object containing the times and statuses of observations for which the metric will be evaluated
#' @param risk ignored, left for compatibility with other metrics
#' @param surv a matrix containing the predicted survival functions for the considered observations, each row represents a single observation, whereas each column one time point
#' @param times a vector of time points at which the survival function was evaluated
#'
#' @return numeric from 0 to 1, lower scores are better (brier score of 0.25 represents a model which returns always returns 0.5 as the predicted survival function)
#'
#' @section References:
#' - \[1\] Brier, Glenn W. ["Verification of forecasts expressed in terms of probability."](https://journals.ametsoc.org/view/journals/mwre/78/1/1520-0493_1950_078_0001_vofeit_2_0_co_2.xml) Monthly Weather Review 78.1 (1950): 1-3.
#' - \[2\] Graf, Erika, et al. ["Assessment and comparison of prognostic classification schemes for survival data."](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0258%2819990915/30%2918%3A17/18%3C2529%3A%3AAID-SIM274%3E3.0.CO%3B2-5) Statistics in Medicine 18.17‐18 (1999): 2529-2545.
#'
#' @rdname brier_score
#' @seealso [cd_auc()]
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' y <- cph_exp$y
#' times <- cph_exp$times
#' surv <- cph_exp$predict_survival_function(cph, cph_exp$data, times)
#'
#' brier_score(y, surv = surv, times = times)
#' loss_brier_score(y, surv = surv, times = times)
#'
#' @export
brier_score <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL) {
    # if times is not provided use
    if (is.null(times)) times <- sort(unique(y_true[, 1]))
    else times <- sort(unique(times))


    # calculate the inverse probability of censoring weights
    y_true[, 2] <- 1 - y_true[, 2]
    km <- survival::survfit(y_true ~ 1)
    G <- stepfun(km$time, c(1, km$surv))

    n_cols <- length(times)
    n_rows <- length(y_true)

    y <- matrix(rep(y_true[, 1], n_cols), ncol = n_cols)
    # return to original deltas!!
    delta <- matrix(rep(1 - y_true[, 2], n_cols), ncol = n_cols)
    ti <- t(matrix(rep(times, n_rows), ncol = n_rows))

    gti <- matrix(G(ti), ncol = n_cols, nrow = n_rows)
    gy  <- matrix(G(y), ncol = n_cols, nrow = n_rows)

    ind_1 <- ifelse(y <= ti & delta == 1, 1, 0)
    ind_2 <- ifelse(y > ti, 1, 0)

    brier_score <- ind_1 * (surv^2) / gy + ind_2 * ((1 - surv)^2) / gti

    apply(brier_score, 2, mean, na.rm = TRUE)

}
attr(brier_score, "loss_name") <- "Brier score"
attr(brier_score, "loss_type") <- "time-dependent"

#' @rdname brier_score
#' @export
loss_brier_score <- brier_score
attr(loss_brier_score, "loss_name") <- "Brier score"
attr(loss_brier_score, "loss_type") <- "time-dependent"

#' Calculate Cumulative/Dynamic AUC
#'
#' This function calculates the Cumulative/Dynamic AUC metric for a survival model. It is done using the
#' estimator proposed proposed by Uno et al. \[[1](https://www.jstor.org/stable/27639883#metadata_info_tab_contents)\],
#' and Hung and Chang \[[2](https://www.jstor.org/stable/41000414#metadata_info_tab_contents)\].
#'
#' C/D AUC is an extension of the AUC metric known from classification models.
#' Its values represent the model's performance at specific time points.
#' It can be integrated over the considered time range.
#'
#' @param y_true a `survival::Surv` object containing the times and statuses of observations for which the metric will be evaluated
#' @param risk ignored, left for compatibility with other metrics
#' @param surv a matrix containing the predicted survival functions for the considered observations, each row represents a single observation, whereas each column one time point
#' @param times a vector of time points at which the survival function was evaluated
#'
#' @return a numeric vector of length equal to the length of the times vector, each value (from the range from 0 to 1) represents the AUC metric at a specific time point, with higher values indicating better performance.
#'
#' @section References:
#' - \[1\] Uno, Hajime, et al. ["Evaluating prediction rules for t-year survivors with censored regression models."](https://www.jstor.org/stable/27639883#metadata_info_tab_contents) Journal of the American Statistical Association 102.478 (2007): 527-537.
#' - \[2\] Hung, Hung, and Chin‐Tsang Chiang. ["Optimal composite markers for time‐dependent receiver operating characteristic curves with censored survival data."](https://www.jstor.org/stable/41000414#metadata_info_tab_contents) Scandinavian Journal of Statistics 37.4 (2010): 664-679.
#'
#' @rdname cd_auc
#' @seealso [loss_one_minus_cd_auc()] [integrated_cd_auc()] [brier_score()]
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' y <- cph_exp$y
#' times <- cph_exp$times
#' surv <- cph_exp$predict_survival_function(cph, cph_exp$data, times)
#'
#' cd_auc(y, surv = surv, times = times)
#'
#' @export
cd_auc <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL) {

    y_true[, 2] <- 1 - y_true[, 2]
    km <- survival::survfit(y_true ~ 1)
    G <- stepfun(km$time, c(1, km$surv))

    n_cols <- length(times)
    n_rows <- length(y_true)

    yi <- matrix(rep(y_true[, 1], n_rows), ncol = n_rows)
    yj <- t(matrix(rep(y_true[, 1], n_rows), ncol = n_rows))

    delta <- matrix(rep(1 - y_true[, 2], n_rows), ncol = n_rows)
    ti <- t(matrix(rep(times, n_rows), ncol = n_rows))

    res <- numeric(length(times))

    results <- lapply(1:length(times), function(tt) {
        time <- times[tt]

        survi <- matrix(rep(surv[, tt], n_rows), ncol = n_rows)
        survj <- t(matrix(rep(surv[, tt], n_rows), ncol = n_rows))

        top <- sum(ifelse(yj > time & yi <= time &  survj > survi, 1, 0) / G(time))

        bl <- sum(ifelse(yi[, 1] > time, 1, 0))
        br <- sum(ifelse(yi[, 1] <= time, 1, 0) / G(time))

        top / (bl * br)
    })

    unlist(results)

}
attr(cd_auc, "loss_name") <- "C/D AUC"
attr(cd_auc, "loss_type") <- "time-dependent"


#' Calculate Cumulative/Dynamic AUC loss
#'
#' This function subtracts the C/D AUC metric from one to obtain a loss function whose lower values indicate better model performance (useful for permutational feature importance)
#'
#' @inheritParams cd_auc
#'
#' @return a numeric vector of length equal to the length of the times vector, each value (from the range from 0 to 1) represents 1 - AUC metric at a specific time point, with lower values indicating better performance.
#'
#' #' @section References:
#' - \[1\] Uno, Hajime, et al. ["Evaluating prediction rules for t-year survivors with censored regression models."](https://www.jstor.org/stable/27639883#metadata_info_tab_contents) Journal of the American Statistical Association 102.478 (2007): 527-537.
#' - \[2\] Hung, Hung, and Chin‐Tsang Chiang. ["Optimal composite markers for time‐dependent receiver operating characteristic curves with censored survival data."](https://www.jstor.org/stable/41000414#metadata_info_tab_contents) Scandinavian Journal of Statistics 37.4 (2010): 664-679.
#'
#' @rdname loss_one_minus_cd_auc
#' @seealso [cd_auc()]
#'
#' @examples
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' y <- cph_exp$y
#' times <- cph_exp$times
#' surv <- cph_exp$predict_survival_function(cph, cph_exp$data, times)
#'
#' loss_one_minus_cd_auc(y, surv = surv, times = times)
#'
#' @export
loss_one_minus_cd_auc <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL) {
    1 - cd_auc(y_true = y_true, risk = risk, surv = surv, times = times)
}
attr(loss_one_minus_cd_auc, "loss_name") <- "One minus C/D AUC"
attr(loss_one_minus_cd_auc, "loss_type") <- "time-dependent"

#' Calculate integrated C/D AUC
#'
#' This function calculates the integrated Cumulative/Dynamic AUC metric for a survival model.
#'
#' It is useful to see how a model performs as a whole, not at specific time points, for example for easier comparison. This function allows for calculating the integral of the C/D AUC metric numerically using the trapezoid method.
#'
#' @param y_true a `survival::Surv` object containing the times and statuses of observations for which the metric will be evaluated
#' @param risk ignored, left for compatibility with other metrics
#' @param surv a matrix containing the predicted survival functions for the considered observations, each row represents a single observation, whereas each column one time point
#' @param times a vector of time points at which the survival function was evaluated
#' @param auc a vector containing already calculated AUC metric at the time points specified in the times parameter. If this is provided all arguments except `times` and `auc` are ignored
#'
#' @return numeric from 0 to 1, higher values indicate better performance
#'
#' #' @section References:
#' - \[1\] Uno, Hajime, et al. ["Evaluating prediction rules for t-year survivors with censored regression models."](https://www.jstor.org/stable/27639883#metadata_info_tab_contents) Journal of the American Statistical Association 102.478 (2007): 527-537.
#' - \[2\] Hung, Hung, and Chin‐Tsang Chiang. ["Optimal composite markers for time‐dependent receiver operating characteristic curves with censored survival data."](https://www.jstor.org/stable/41000414#metadata_info_tab_contents) Scandinavian Journal of Statistics 37.4 (2010): 664-679.
#'
#' @rdname integrated_cd_auc
#' @seealso [cd_auc()] [loss_one_minus_cd_auc()]
#'
#' @examples
#'
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' y <- cph_exp$y
#' times <- cph_exp$times
#' surv <- cph_exp$predict_survival_function(cph, cph_exp$data, times)
#'
#'
#' # calculating directly
#' integrated_cd_auc(y, surv = surv, times = times)
#'
#' # calculating based on given auc vector
#' auc <- cd_auc(y, surv = surv, times = times)
#' integrated_cd_auc(times = times, auc = auc)
#'
#' @export
integrated_cd_auc <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL, auc = NULL) {
    if (is.null(auc)) {
        auc <- cd_auc(y_true, risk, surv, times)
    }

    auc <- as.numeric(auc)
    ind_to_drop <- is.na(auc)
    times <- times[!ind_to_drop]
    auc <- auc[!ind_to_drop]
    n <- length(auc)

    iauc <- (auc[1:(n - 1)] + auc[2:n]) * diff(times) / 2

    cumsum(c(0, iauc))[length(cumsum(c(0, iauc)))] / (max(times) - min(times))
}
attr(integrated_cd_auc, "loss_name") <- "integrated C/D AUC"
attr(integrated_cd_auc, "loss_type") <- "integrated"



#' Calculate integrated C/D AUC loss
#'
#' This function subtracts integrated the C/D AUC metric from one to obtain a loss function whose lower values indicate better model performance (useful for permutational feature importance)
#'
#' @inheritParams integrated_cd_auc
#'
#' @return numeric from 0 to 1, lower values indicate better performance
#'
#' #' @section References:
#' - \[1\] Uno, Hajime, et al. ["Evaluating prediction rules for t-year survivors with censored regression models."](https://www.jstor.org/stable/27639883#metadata_info_tab_contents) Journal of the American Statistical Association 102.478 (2007): 527-537.
#' - \[2\] Hung, Hung, and Chin‐Tsang Chiang. ["Optimal composite markers for time‐dependent receiver operating characteristic curves with censored survival data."](https://www.jstor.org/stable/41000414#metadata_info_tab_contents) Scandinavian Journal of Statistics 37.4 (2010): 664-679.
#'
#' @rdname loss_one_minus_integrated_cd_auc
#' @seealso [integrated_cd_auc()] [cd_auc()] [loss_one_minus_cd_auc()]
#'
#' @examples
#'
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' y <- cph_exp$y
#' times <- cph_exp$times
#' surv <- cph_exp$predict_survival_function(cph, cph_exp$data, times)
#'
#'
#' # calculating directly
#' loss_one_minus_integrated_cd_auc(y, surv = surv, times = times)
#'
#' # calculating based on given auc vector
#' auc <- cd_auc(y, surv = surv, times = times)
#' loss_one_minus_integrated_cd_auc(times = times, auc = auc)
#'
#' @export
loss_one_minus_integrated_cd_auc <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL, auc = NULL) {
    1 - integrated_cd_auc(y_true = y_true, risk = risk, surv = surv, times = times, auc = auc)
}
attr(loss_one_minus_integrated_cd_auc, "loss_name") <- "One minus integrated C/D AUC"
attr(loss_one_minus_integrated_cd_auc, "loss_type") <- "integrated"



#' Calculate integrated Brier score
#'
#' This function calculates the integrated Brier score metric for a survival model.
#'
#' It is useful to see how a model performs as a whole, not at specific time points, for example for easier comparison. This function allows for calculating the integral of Brier score metric numerically using the trapezoid method.
#'
#' @param y_true a `survival::Surv` object containing the times and statuses of observations for which the metric will be evaluated
#' @param risk ignored, left for compatibility with other metrics
#' @param surv a matrix containing the predicted survival functions for the considered observations, each row represents a single observation, whereas each column one time point
#' @param times a vector of time points at which the survival function was evaluated
#' @param brier a vector containing already calculated Brier score metric at the time points specified in the times parameter. If this is provided all arguments except `times` and `brier` are ignored
#'
#' @return numeric from 0 to 1, lower values indicate better performance
#'
#' @section References:
#' - \[1\] Brier, Glenn W. ["Verification of forecasts expressed in terms of probability."](https://journals.ametsoc.org/view/journals/mwre/78/1/1520-0493_1950_078_0001_vofeit_2_0_co_2.xml) Monthly Weather Review 78.1 (1950): 1-3.
#' - \[2\] Graf, Erika, et al. ["Assessment and comparison of prognostic classification schemes for survival data."](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0258%2819990915/30%2918%3A17/18%3C2529%3A%3AAID-SIM274%3E3.0.CO%3B2-5) Statistics in Medicine 18.17‐18 (1999): 2529-2545.
#'
#' @rdname integrated_brier_score
#' @seealso [brier_score()] [integrated_cd_auc()] [loss_one_minus_integrated_cd_auc()]
#'
#' @examples
#'
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' y <- cph_exp$y
#' times <- cph_exp$times
#' surv <- cph_exp$predict_survival_function(cph, cph_exp$data, times)
#'
#'
#' # calculating directly
#' integrated_brier_score(y, surv = surv, times = times)
#'
#' # calculating based on given auc vector
#' brier_score <- brier_score(y, surv = surv, times = times)
#' integrated_brier_score(times = times, brier = brier_score)
#'
#' @export
integrated_brier_score <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL, brier = NULL) {

    if (is.null(brier)) {
        brier <- brier_score(y_true, risk, surv, times)
    }

    brier <- as.numeric(brier)
    ind_to_drop <- is.na(brier)
    times <- times[!ind_to_drop]
    brier <- brier[!ind_to_drop]
    n <- length(brier)

    ibs <- (brier[1:(n - 1)] + brier[2:n]) * diff(times) / 2

    cumsum(c(0, ibs))[length(cumsum(c(0, ibs)))] / (max(times) - min(times))

}
attr(integrated_brier_score, "loss_name") <- "integrated Brier score"

#' @rdname integrated_brier_score
#' @export
loss_integrated_brier_score <- integrated_brier_score
attr(loss_integrated_brier_score, "loss_name") <- "integrated Brier score"
attr(loss_integrated_brier_score, "loss_type") <- "integrated"

#' Adapt mlr3proba measures for use with survex
#'
#' This function allows for usage of standardized measures from the mlr3proba package with `survex`.
#'
#' @param measure - a `MeasureSurv` object from the `mlr3proba` package, the object to adapt
#' @param reverse - boolean, FALSE by default, whether the metric should be reversed in order to be treated as loss (for permutational variable importance we need functions with lower values indicating better performance). If TRUE, the new metric value will be (1 - metric_value)
#' @param ... - other parameters, currently ignored
#'
#' @return a function with standardized parameters (`y_true`, `risk`, `surv`, `times`) that can be used to calculate loss
#'
#' if(FALSE){
#'   measure <- msr("surv.calib_beta")
#'   mlr_measure <- loss_adapt_mlr3proba(measure)
#' }
#'
#' @export
loss_adapt_mlr3proba <- function(measure, reverse = FALSE, ...){

    loss_function <- function(y_true = NULL, risk = NULL, surv = NULL, times = NULL){

        colnames(surv) <- times

        surv_pred <- PredictionSurv$new(
            row_ids = 1:length(y_true),
            truth = y_true,
            crank = risk,
            distr = surv,
            task = list(truth = y_true)
        )

        output <- surv_pred$score(measure)
        names(output) <- NULL

        if (reverse) output <- (1 - output)

        return(output)
    }

    if (reverse) attr(loss_function, "loss_name") <- paste("one minus", measure$id)
    else attr(loss_function, "loss_name") <- measure$id
    attr(loss_function, "loss_type") <- "integrated"

    return(loss_function)
}




