#' Instance Level Parts of Survival Model Predictions
#'
#' This function decomposes the model prediction into individual parts, which are attributions of particular variables. The explanations can be made via the SurvLIME and SurvSHAP(t) methods.
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param new_observation a new observation for which prediction need to be explained
#' @param ... other parameters which are passed to `iBreakDown::break_down` if `output_type=="risk"`, or if `output_type=="survival"` to `surv_shap()` or `surv_lime()` functions depending on the selected type
#' @param N the maximum number of observations used for calculation of attributions. If `NULL` (default) all observations will be used.
#' @param type if `output_type == "survival"` must be either `"survshap"` or `"survlime"`, otherwise refer to the `DALEX::predict_parts`
#' @param output_type either `"survival"` or `"risk"` the type of survival model output that should be considered for explanations. If `"survival"` the explanations are based on the survival function. Otherwise the scalar risk predictions are used by the `DALEX::predict_parts` function.
#'
#' @return An object of class `"predict_parts_survival"` and additional classes depending on the type of explanations. It is a list with the element `result` containing the results of the calculation.
#'
#'
#' @section Additional parameters:
#'
#' There are additional parameters that are passed to internal functions
#'
#' * for `surv_lime`
#'     * `N` -  a positive integer, number of observations generated in the neighbourhood
#'     * `distance_metric` -  character, name of the distance metric to be used, only `"euclidean"` is implemented
#'     * `kernel_width` -  a numeric, parameter used for calculating weights, by default it's `sqrt(ncol(data)*0.75)`
#'     * `sampling_method` -  character, name of the method of generating neighbourhood, only `"gaussian"` is implemented
#'     * `sample_around_instance` -  logical, if the neighbourhood should be generated with the new observation as the center (default), or should the mean of the whole dataset be used as the center
#'     * `max_iter` -  a numeric, maximal number of iteration for the optimization problem
#'     * `categorical_variables` -  character vector, names of variables that should be treated as categories (factors are included by default)
#'     * `k` -  a small positive number > 1, added to chf before taking log, so that weigths aren't negative
#' * for `surv_shap`
#'     * `timestamps` -  a numeric vector, time points at which the survival function will be evaluated
#'     * `y_true` -  a two element numeric vector or matrix of one row and two columns, the first element being the true observed time and the second the status of the observation, used for plotting
#'     * `calculation_method` -  a character, only `"kernel"` is implemented for now.
#'     * `aggregation_method` -  a character, either `"mean_absolute"` or `"integral"`, `"max_absolute"`, `"sum_of_squares"`
#'
#' @section References:
#' - [Krzyziński, Mateusz, et al. "SurvSHAP (t): Time-dependent explanations of machine learning survival models." arXiv preprint arXiv:2208.11080 (2022).](https://arxiv.org/abs/2208.11080)
#' - [Kovalev, Maxim S., Lev V. Utkin, and Ernest M. Kasimov. "SurvLIME: A method for explaining machine learning survival models." Knowledge-Based Systems 203 (2020): 106164.](https://www.sciencedirect.com/science/article/pii/S0950705120304044?casa_token=6e9cyk_ji3AAAAAA:tbqo33MsZvNC9nrSGabZdLfPtZTsvsvZTHYQCM2aEhumLI5D46U7ovhr37EaYUhmKZrw45JzDhg)
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' cph_exp <- explain(cph)
#'
#' cph_predict_parts_survshap <- predict_parts(cph_exp, new_observation = veteran[1, -c(3, 4)])
#' head(cph_predict_parts_survshap$result)
#' plot(cph_predict_parts_survshap)
#'
#' cph_predict_parts_survlime <- predict_parts(cph_exp, new_observation = veteran[1, -c(3, 4)],
#'                                             type = "survlime")
#' head(cph_predict_parts_survlime$result)
#' plot(cph_predict_parts_survlime, type = "local_importance")
#' }
#'
#' @rdname predict_parts.surv_explainer
#' @export
predict_parts <- function(explainer, ...) UseMethod("predict_parts", explainer)

#' @rdname predict_parts.surv_explainer
#' @export
predict_parts.surv_explainer <- function(explainer, new_observation, ..., N = NULL, type = "survshap", output_type = "survival") {


    label <- explainer$label

    if (output_type == "risk") {
        return (DALEX::predict_parts(explainer = explainer,
                             new_observation = new_observation,
                             ... = ...,
                             N = N,
                             type = type))
    }
    else {

        res <- switch(type,
               "survshap" = surv_shap(explainer, new_observation, ...),
               "survlime" = surv_lime(explainer, new_observation, ...),
               stop("Only `survshap` and `survlime` methods are implemented for now"))

    }

    attr(res, "label") <- label
    class(res) <- c('predict_parts_survival', class(res))
    res

}


#' @export
predict_parts.default <- DALEX::predict_parts
