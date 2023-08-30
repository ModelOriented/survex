#' Dataset Level Model Diagnostics
#'
#' This function calculates martingale and deviance residuals.
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#'
#' @return An object of class `c("model_diagnostics_survival")`. It's a list with the explanations in the `result` element.
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
#' cph_residuals <- model_diagnostics(cph_exp)
#' rsf_residuals <- model_diagnostics(rsf_ranger_exp)
#'
#' head(cph_residuals$result)
#' plot(cph_residuals, rsf_residuals, xvariable = "age")
#' plot(cph_residuals, rsf_residuals, plot_type = "Cox-Snell")
#'
#' }
#' @rdname model_diagnostics.surv_explainer
#' @export
model_diagnostics <- function(explainer) UseMethod("model_diagnostics", explainer)

#' @rdname model_diagnostics.surv_explainer
#' @export
model_diagnostics.surv_explainer <- function(explainer) {
    test_explainer(explainer,  has_data = TRUE, has_y = TRUE, has_chf = TRUE, function_name = "model_diagnostics")

    n <- nrow(explainer$data)
    original_times <- explainer$y[, 1]
    statuses <-  explainer$y[, 2]

    unique_times <- sort(unique(original_times))
    which_el <- matrix(c(1:n, match(original_times, unique_times)),
                      nrow = n)
    chf_preds <-
       predict(explainer, times = unique_times, output_type = "chf")
    cox_snell_residuals <- chf_preds[which_el]
    martingale_residuals <- statuses - cox_snell_residuals
    deviance_residuals <- sign(martingale_residuals) *
       sqrt(-2 * (
           martingale_residuals + statuses *
               log(statuses - martingale_residuals)
       ))
    # cox_snell_residuals[statuses == 0] <-
    #     cox_snell_residuals[statuses == 0] + 1
    # modification for censored observations

    result <- cbind(
       data.frame(
           "time" = original_times,
           "status" = factor(statuses),
           "cox_snell_residuals" = cox_snell_residuals,
           "martingale_residuals" = martingale_residuals,
           "deviance_residuals" = deviance_residuals,
           "label" = explainer$label
       ),
       explainer$data
    )

    res <- list(result = result)
    class(res) <- c("model_diagnostics_survival", class(res))
    res
}
