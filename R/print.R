#' @export
print.model_profile_survival <- function(x, ...) {
    res <- x$result

    res <- res[order(res$`_vname_`), ]
    text <- paste0("Partial dependence for the ", unique(res$`_label_`), " model:\n")
    cat(text)
    print.data.frame(res[, !colnames(res) %in% c("_label_", "_ids_")], ...)
}


#' @export
print.surv_ceteris_paribus <- function(x, ...) {
    res <- x$result
    cat("Ceteris paribus for observation:\n\n")
    print.data.frame(x$variable_values, row.names = FALSE)
    cat("\n")
    print(res, ...)
}


#' @export
print.surv_feature_importance <- function(x, ...) {

    res <- x$result
    text <- paste0("Permutational feature importance for the ", unique(res$label), " model:\n")
    cat(text)
    print.data.frame(x$result[res$`_permutation_` == 0, !colnames(res) %in% c("_baseline_", "_permutation_", "label")], ...)
}


#' @export
print.surv_lime <- function(x, ...) {
    res <- x$result

    print_result <- rbind(beta = res, `X` = x$variable_values, `local importance (X*beta)` = res * as.numeric(x$variable_values))
    cat("SurvLIME for observation:\n\n")
    print.data.frame(x$variable_values, row.names = FALSE)
    cat("\n")
    print(t(print_result))
}


#' @export
print.surv_shap <- function(x, ...) {

    res <- x$result
    cat("SurvSHAP(t) for observation:\n\n")
    print.data.frame(x$variable_values, row.names = FALSE)
    cat("\n")
    print(res, ...)
}
