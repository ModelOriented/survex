#' Helper functions for `model_profile.R`
#'
#' @param x an object containing calculated ceteris_paribus profiles
#' @param ... other parameters, ignored
#' @param variable_type character, either `"numerical"` or `"categorical"`, the type of variable to be calculated, if left `NULL` (default), both are calculated
#' @param groups unused, left for compatibility
#' @param type character, only `"partial"` is implemented
#' @param variables a character vector containing names of variables to be explained
#' @param center logical, if the profiles should be centered before aggregations
#'
#' @return A data.frame with calculated results.
#'
#' @keywords internal
#' @importFrom stats na.omit aggregate
surv_aggregate_profiles <- function(x, ...,
                                         variable_type = NULL,
                                         groups = NULL,
                                         type = "partial",
                                         variables = NULL,
                                         center = FALSE) {




    all_profiles <- x$result
    class(all_profiles) <- "data.frame"

    all_profiles$`_ids_` <- factor(all_profiles$`_ids_`)


    # variables to use
    all_variables <- na.omit(as.character(unique(all_profiles$`_vname_`)))
    if (!is.null(variables)) {
        all_variables_intersect <- intersect(all_variables, variables)
        if (length(all_variables_intersect) == 0) stop(paste0("parameter variables do not overlap with ", paste(all_variables, collapse = ", ")))
        all_variables <- all_variables_intersect
    }

    if (!is.null(variable_type) && variable_type == "numerical") {
        all_profiles <- all_profiles[all_profiles$`_vtype_` == "numerical", ]
    }

    if (!is.null(variable_type) && variable_type == "categorical") {
        all_profiles <- all_profiles[all_profiles$`_vtype_` == "categorical", ]
    }

    all_variables <- intersect(all_variables, unique(all_profiles$`_vname_`))

    # select only suitable variables
    all_profiles <- all_profiles[all_profiles$`_vname_` %in% all_variables, ]
    # create _x_
    tmp <- as.character(all_profiles$`_vname_`)
    for (viname in unique(tmp)) {
        all_profiles$`_x_`[tmp == viname] <- all_profiles[tmp == viname, viname]
    }

    if (!inherits(class(all_profiles), "data.frame")) {
        all_profiles <- as.data.frame(all_profiles)
    }

    # change x column to proper character values
    for (variable in all_variables) {
        if (variable %in% all_profiles[all_profiles$`_vtype_` == "categorical", "_vname_"])
        all_profiles[all_profiles$`_vname_` == variable, ]$`_x_` <- as.character(apply(all_profiles[all_profiles$`_vname_` == variable, ], 1, function(all_profiles) all_profiles[all_profiles["_vname_"]]))
    }

    if (type == "partial") {
        aggregated_profiles <- surv_aggregate_profiles_partial(all_profiles)
        class(aggregated_profiles) <- c("aggregated_survival_profiles_explainer",
                                        "partial_dependence_survival_explainer",
                                        "data.frame")
    }
    if (type == "conditional") {
        stop("Not implemented")
    }
    if (type == "accumulated") {
        stop("Not implemented")
    }

    return(aggregated_profiles)
}


surv_aggregate_profiles_partial <- function(all_profiles) {

    tmp <- all_profiles[, c("_vname_", "_vtype_", "_label_", "_x_", "_yhat_", "_times_")]
    aggregated_profiles <- aggregate(tmp$`_yhat_`, by = list(tmp$`_vname_`, tmp$`_vtype_`, tmp$`_label_`, tmp$`_x_`, tmp$`_times_`), FUN = mean, na.rm = TRUE)
    colnames(aggregated_profiles) <- c("_vname_",  "_vtype_", "_label_", "_x_", "_times_", "_yhat_")
    aggregated_profiles$`_ids_` <- 0

    # for factors, keep proper order
    # as in https://github.com/ModelOriented/ingredients/issues/82
    if (!is.numeric(all_profiles$`_x_`)) {
        aggregated_profiles$`_x_` <- factor(aggregated_profiles$`_x_`, levels = unique(all_profiles$`_x_`))
        aggregated_profiles <- aggregated_profiles[order(aggregated_profiles$`_x_`), ]
    }

    aggregated_profiles


}
