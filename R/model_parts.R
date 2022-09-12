#' Dataset Level Variable Importance for Survival Models
#'
#' This function calculates variable importance as a change in the loss function after the variable values permutations.
#'
#'
#' @param explainer a model to be explained, preprocessed by the `explain()` function
#' @param loss_function a function that will be used to assess variable importance, by default `loss_brier_score` for survival models. The function can be supplied manually but has to have these named parameters (`y_true`, `risk`, `surv`, `times`), where `y_true` represents the `survival::Surv` object with observed times and statuses, `risk` is the risk score calculated by the model, and `surv` is the survival function for each observation evaluated at `times`.
#' @param ... other parameters passed to `DALEX::model_parts` if `output_type == "risk"`, otherwise passed to internal functions.
#' @param type a character vector, if `"raw"` the results are losses after the permutation, if `"ratio"` the results are in the form `loss/loss_full_model` and if `"difference"` the results are of the form `loss - loss_full_model`
#' @param output_type either `"survival"` or `"risk"` the type of survival model output that should be used for explanations. If `"survival"` the explanations are based on the survival function. Otherwise the scalar risk predictions are used by the `DALEX::model_profile` function.
#' @param N number of observations that should be sampled for calculation of variable importance. If `NULL` then variable importance will be calculated on the whole dataset.
#' @inheritDotParams surv_feature_importance -x
#' @inheritDotParams surv_integrated_feature_importance -x
#'
#'
#' @return An object of class `c("model_parts_survival", "surv_feature_importance")`. It's a list with the explanations in the `result` element.
#'
#' @details
#' *Note*: This function can be run within `progressr::with_progress()` to display a progress bar, as the execution can take long, especially on large datasets.
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' rsf_ranger <- ranger::ranger(Surv(time, status) ~ .,
#'   data = veteran,
#'   respect.unordered.factors = TRUE,
#'   num.trees = 100,
#'   mtry = 3,
#'   max.depth = 5
#' )
#'
#' cph_exp <- explain(cph)
#'
#' rsf_ranger_exp <- explain(rsf_ranger,
#'   data = veteran[, -c(3, 4)],
#'   y = Surv(veteran$time, veteran$status)
#' )
#'
#' cph_model_parts_brier <- model_parts(cph_exp)
#' print(head(cph_model_parts_brier$result))
#' plot(cph_model_parts_brier)
#'
#'
#' rsf_ranger_model_parts <- model_parts(rsf_ranger_exp)
#' print(head(rsf_ranger_model_parts$result))
#' plot(cph_model_parts_brier, rsf_ranger_model_parts)
#' }
#' @rdname model_parts.surv_explainer
#' @export
model_parts <- function(explainer, ...) UseMethod("model_parts", explainer)

#' @export
model_parts.surv_explainer <- function(explainer,
                                       loss_function = survex::loss_brier_score,
                                       ...,
                                       type = "raw",
                                       output_type = "survival",
                                       N = 1000) {
  if (!(type %in% c("difference", "ratio", "raw", "variable_importance"))) stop("Type shall be one of 'variable_importance', 'difference', 'ratio', 'raw'")
  if (!(output_type %in% c("survival", "risk"))) stop("output_type should be 'survival' or 'risk'")
  if (type == "variable_importance") type <- "raw" # it's an alias

  switch(output_type,
    "risk" = DALEX::model_parts(
      explainer = explainer,
      loss_function = loss_function,
      ... = ...,
      type = type,
      N = N
    ),
    "survival" = {
      test_explainer(explainer, has_data = TRUE, has_y = TRUE, has_survival = TRUE, function_name = "model_parts")

      if (attr(loss_function, "loss_name") %in% c("integrated Brier score", "One minus integrated C/D AUC", "One minus C-Index")) {
        res <- surv_integrated_feature_importance(
          x = explainer,
          loss_function = loss_function,
          type = type,
          N = N,
          ...
        )
        class(res) <- c("model_parts", class(res))
        return(res)
      } else {
        res <- surv_feature_importance(
          x = explainer,
          loss_function = loss_function,
          type = type,
          N = N,
          ...
        )
        class(res) <- c("model_parts_survival", class(res))
        res
      }
    },
    stop("Type should be either `survival` or `risk`")
  )
}

#' @export
model_parts.default <- DALEX::model_parts
