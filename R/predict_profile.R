#' Instance Level Profile as Ceteris Paribus for Survival Models
#'
#' This function calculates Ceteris Paribus Proifles for a specific observation with the possibility to take the time dimension into account.
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param new_observation a new observation for which the prediction need to be explained
#' @param variables a character vector containing names of variables to be explained
#' @param categorical_variables a character vector of names of additional variables which should be treated as categorical (factors are automatically treated as categorical variables). If it contains variable names not present in the `variables` argument, they will be added at the end.
#' @param ... additional parameters passed to `DALEX::predict_profile` if `output_type =="risk"`
#' @param type character, only `"ceteris_paribus"` is implemented
#' @param output_type either `"survival"` or `"risk"` the type of survival model output that should be considered for explanations. If `"survival"` the explanations are based on the survival function. Otherwise the scalar risk predictions are used by the `DALEX::predict_profile` function.
#' @param variable_splits_type character, decides how variable grids should be calculated. Use `"quantiles"` for percentiles or `"uniform"` (default) to get uniform grid of points.
#'
#' @return An object of class `c("predict_profile_survival", "surv_ceteris_paribus")`. It is a list with the final result in the `result` element.
#'
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
#' rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#'
#' cph_exp <- explain(cph)
#' rsf_src_exp <- explain(rsf_src)
#'
#' cph_predict_profile <- predict_profile(cph_exp, veteran[2, -c(3, 4)],
#'                                        variables = c("trt", "celltype", "karno", "age"),
#'                                        categorical_variables = "trt")
#' plot(cph_predict_profile, facet_ncol = 2)
#'
#'
#' rsf_predict_profile <- predict_profile(rsf_src_exp, veteran[5, -c(3, 4)], variables = "karno")
#' plot(cph_predict_profile, numerical_plot_type = "contours")
#' }
#'
#' @rdname predict_profile.surv_explainer
#' @export
predict_profile <- function(explainer,
                            new_observation,
                            variables = NULL,
                            categorical_variables = NULL,
                            ...,
                            type = "ceteris_paribus",
                            variable_splits_type = "uniform")
    UseMethod("predict_profile", explainer)

#' @rdname predict_profile.surv_explainer
#' @export
predict_profile.surv_explainer <- function(explainer,
                                           new_observation,
                                           variables = NULL,
                                           categorical_variables = NULL,
                                           ...,
                                           type = "ceteris_paribus",
                                           output_type = "survival",
                                           variable_splits_type = "uniform")
{

    variables <- unique(variables, categorical_variables)
    if (!type %in% "ceteris_paribus") stop("Type not supported")
    if (!output_type %in% c("risk", "survival")) stop("output_type not supported")

    if (output_type == "risk") {
        return(DALEX::predict_profile(explainer = explainer,
                                      new_observation = new_observation,
                                      variables = variables,
                                      ... = ...,
                                      type = type,
                                      variable_splits_type = variable_splits_type))

    }
    if (output_type == "survival") {
        if (type == "ceteris_paribus") {
            res <- surv_ceteris_paribus(explainer,
                                        new_observation = new_observation,
                                        variables = variables,
                                        categorical_variables = categorical_variables,
                                        variable_splits_type = variable_splits_type,
                                        ...)
            class(res) <- c("predict_profile_survival", class(res))
            return(res)
        }
        else stop("For survival output only type=`ceteris_paribus` is implemented")
    }

}

#' @export
predict_profile.default <- DALEX::predict_profile
