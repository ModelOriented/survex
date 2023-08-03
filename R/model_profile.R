#' Dataset Level Variable Profile as Partial Dependence Explanations for Survival Models
#'
#' This function calculates explanations on a dataset level that help explore model response as a function of selected variables.
#' The explanations are calculated as an extension of Partial Dependence Profiles with the inclusion of the time dimension.
#'
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param variables character, a vector of names of variables to be explained
#' @param N number of observations used for the calculation of aggregated profiles. By default `100`. If `NULL` all observations are used.
#' @param ... other parameters passed to `DALEX::model_profile` if `output_type == "risk"`, otherwise ignored
#' @param categorical_variables character, a vector of names of additional variables which should be treated as categorical (factors are automatically treated as categorical variables). If it contains variable names not present in the `variables` argument, they will be added at the end.
#' @param grid_points maximum number of points for profile calculations. Note that the final number of points may be lower than grid_points. Will be passed to internal function. By default `51`.
#' @param variable_splits_type character, decides how variable grids should be calculated. Use `"quantiles"` for percentiles or `"uniform"` (default) to get uniform grid of points.
#' @param groups if `output_type == "risk"` a variable name that will be used for grouping. By default `NULL`, so no groups are calculated. If `output_type == "survival"` then ignored
#' @param k passed to `DALEX::model_profile` if `output_type == "risk"`, otherwise ignored
#' @param center logical, should profiles be centered before clustering
#' @param type the type of variable profile, `"partial"` for Partial Dependence, `"accumulated"` for Accumulated Local Effects, or `"conditional"` (available only for `output_type == "risk"`)
#' @param output_type either `"survival"` or `"risk"` the type of survival model output that should be considered for explanations. If `"survival"` the explanations are based on the survival function. Otherwise the scalar risk predictions are used by the `DALEX::model_profile` function.
#'
#' @return An object of class `model_profile_survival`. It is a list with the element `result` containing the results of the calculation.
#'
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
#' cph_model_profile <- model_profile(cph_exp, output_type = "survival",
#'                                    variables = c("age"))
#'
#' head(cph_model_profile$result)
#'
#' plot(cph_model_profile)
#'
#' rsf_model_profile <- model_profile(rsf_src_exp, output_type = "survival",
#'                                    variables = c("age", "celltype"))
#'
#' head(rsf_model_profile$result)
#'
#' plot(rsf_model_profile, variables = c("age", "celltype"), numerical_plot_type = "contours")
#' }
#'
#' @rdname model_profile.surv_explainer
#' @export
model_profile <- function(explainer,
                          variables = NULL,
                          N = 100,
                          ...,
                          groups = NULL,
                          k = NULL,
                          center = FALSE,
                          type = "partial",
                          output_type = "survival") UseMethod("model_profile", explainer)

#' @rdname model_profile.surv_explainer
#' @export
model_profile.surv_explainer <- function(explainer,
                                         variables = NULL,
                                         N = 100,
                                         ...,
                                         categorical_variables = NULL,
                                         grid_points = 51,
                                         variable_splits_type = "uniform",
                                         groups = NULL,
                                         k = NULL,
                                         center = TRUE,
                                         type = "partial",
                                         output_type = "survival") {

    variables <- unique(variables, categorical_variables)
    switch(output_type,
            "risk" = DALEX::model_profile(explainer = explainer,
                                                   variables = variables,
                                                   N = N,
                                                   ... = ...,
                                                   groups = groups,
                                                   k = k,
                                                   center = center,
                                                   type = type),
            "survival" = {
                test_explainer(explainer, "model_profile", has_data = TRUE, has_survival = TRUE)

                data <- explainer$data
                if (!is.null(N) && N < nrow(data)) {
                    ndata <- data[sample(1:nrow(data), N), ]
                } else {
                    ndata <- data
                }

                if (type == "partial"){
                    cp_profiles <- surv_ceteris_paribus(explainer,
                                                        new_observation = ndata,
                                                        variables = variables,
                                                        categorical_variables = categorical_variables,
                                                        grid_points = grid_points,
                                                        variable_splits_type = variable_splits_type,
                                                        ...)

                    result <- surv_aggregate_profiles(cp_profiles, ...,
                                                            variables = variables,
                                                            center = center)
                } else if (type == "accumulated"){
                    cp_profiles <- NULL

                    result <- surv_ale(explainer,
                                       data = ndata,
                                       variables = variables,
                                       categorical_variables = categorical_variables,
                                       grid_points = grid_points,
                                       ...)
                } else {
                    stop("Currently only `partial` and `accumulated` types are implemented")
                }

                ret <- list(eval_times = unique(result$`_times_`),
                            cp_profiles = cp_profiles,
                            result = result,
                            type = type)
                class(ret) <- c("model_profile_survival", "list")
                ret$event_times <- explainer$y[explainer$y[, 1] <= max(explainer$times), 1]
                ret$event_statuses <- explainer$y[explainer$y[, 1] <= max(explainer$times), 2]
                ret
                },

        stop("Currently only `risk` and `survival` output types are implemented")
    )
}

#' @export
model_profile.default <- DALEX::model_profile
