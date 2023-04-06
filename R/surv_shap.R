#' Helper functions for `predict_parts.R`
#'
#' @param explainer an explainer object - model preprocessed by the `explain()` function
#' @param new_observation a new observation for which predictions need to be explained
#' @param ... additional parameters, passed to internal functions
#' @param y_true a two element numeric vector or matrix of one row and two columns, the first element being the true observed time and the second the status of the observation, used for plotting
#' @param calculation_method a character, either `"kernelshap"` for use of `kernelshap` library (providing faster Kernel SHAP with refinements), `"exact_kernel"` for exact Kernel SHAP estimation,
#'   or `"treeshap"` for use of `treeshap` library (efficient implementation to compute SHAP values for tree-based models).
#' @param aggregation_method a character, either `"mean_absolute"` or `"integral"` (default), `"max_absolute"`, `"sum_of_squares"`
#'
#' @return A list, containing the calculated SurvSHAP(t) results in the `result` field
#'
#' @section References:
#' - \[1\] Krzyziński, Mateusz, et al. ["SurvSHAP(t): Time-dependent explanations of machine learning survival models."](https://www.sciencedirect.com/science/article/pii/S0950705122013302) Knowledge-Based Systems 262 (2023): 110234
#'
#' @keywords internal
surv_shap <- function(explainer,
                      new_observation,
                      ...,
                      y_true = NULL,

                      calculation_method = c("kernelshap", "exact_kernel", "treeshap"),
                      aggregation_method = c("integral", "mean_absolute", "max_absolute", "sum_of_squares"),
                      path = "average",
                      B = 25,
                      exact = FALSE
) {
    calculation_method <- match.arg(calculation_method)
    aggregation_method <- match.arg(aggregation_method)

    # make this code work for multiple observations
    stopifnot(ifelse(!is.null(y_true),
                     ifelse(is.matrix(y_true),
                            nrow(new_observation) == nrow(y_true),
                            is.null(dim(y_true)) && length(y_true) == 2L),
                     TRUE))

    if (calculation_method == "kernelshap") {
        if (!requireNamespace("kernelshap", quietly = TRUE)) {
            stop(
                paste0(
                    "Package \"kernelshap\" must be installed to use ",
                    "'calculation_method = \"kernelshap\"'."
                ),
                call. = FALSE
            )
        }
    }
    if (calculation_method == "treeshap") {
        if (!requireNamespace("treeshap", quietly = TRUE)) {
            stop(
                paste0(
                    "Package \"treeshap\" must be installed to use ",
                    "'calculation_method = \"treeshap\"'."
                ),
                call. = FALSE
            )
        }
    }

    test_explainer(explainer, "surv_shap", has_data = TRUE, has_y = TRUE, has_survival = TRUE)

    # make this code also work for 1-row matrix
    col_index <- which(colnames(new_observation) %in% colnames(explainer$data))
    if (is.matrix(new_observation) && nrow(new_observation) == 1) {
        new_observation <- as.matrix(t(new_observation[, col_index]))
    } else {
        new_observation <- new_observation[, col_index]
    }

    if (ncol(explainer$data) != ncol(new_observation)) stop("New observation and data have different number of columns (variables)")

    if (!is.null(y_true)) {
        if (is.matrix(y_true)) {
            # above, we have already checked that nrows of observations are
            # identical to nrows of y_true; thus we do not need to index
            # the first row here
            y_true_ind <- y_true[, 2]
            y_true_time <- y_true[, 1]
        } else {
            y_true_ind <- y_true[2]
            y_true_time <- y_true[1]
        }
    }

    # hack to use rf-model's death times as explainer death times, as
    # treeshap::ranger_surv_fun.unify extracts survival time-points directly
    # from the ranger object for calculating the predictions
    if (calculation_method == "treeshap") {
        if (inherits(explainer$model, "ranger")) {
            explainer$times <- explainer$model$unique.death.times
        } else {
            stop("Calculation method `treeshap` is currently only implemented for `ranger`.")
        }
    }

    res <- list()
    res$eval_times <- explainer$times
    # to display final object correctly, when is.matrix(new_observation) == TRUE
    res$variable_values <- as.data.frame(new_observation)

    res$result <- switch(calculation_method,
                         "exact_kernel" = shap_kernel(explainer, new_observation, ...),
                         "kernelshap" = use_kernelshap(explainer, new_observation, ...),
                         "treeshap" = use_treeshap(explainer, new_observation, ...),
                         stop("Only `exact_kernel`, `kernelshap` and `treeshap` calculation methods are implemented"))

    if (!is.null(y_true)) res$y_true <- c(y_true_time = y_true_time, y_true_ind = y_true_ind)

    res$aggregate <- aggregate_surv_shap(res, aggregation_method)

    class(res) <- "surv_shap"
    res
}


shap_kernel <- function(explainer, new_observation, ...) {

    timestamps <- explainer$times
    p <- ncol(explainer$data)

    target_sf <- explainer$predict_survival_function(explainer$model, new_observation, timestamps)
    sfs <- explainer$predict_survival_function(explainer$model, explainer$data, timestamps)
    baseline_sf <- apply(sfs, 2, mean)


    permutations <- expand.grid(rep(list(0:1), p))
    kernel_weights <- generate_shap_kernel_weights(permutations, p)

    shap_values <- calculate_shap_values(explainer, explainer$model, baseline_sf, explainer$data, permutations, kernel_weights, new_observation, timestamps)

    shap_values <- as.data.frame(shap_values, row.names = colnames(explainer$data))
    colnames(shap_values) <- paste("t=", timestamps, sep = "")
    return (t(shap_values))
}

generate_shap_kernel_weights <- function(permutations, p) {

    apply(permutations, 1, function(row) {
        row <- as.numeric(row)
        num_available_variables = sum(row != 0)

        if (num_available_variables == 0 || num_available_variables == p) 1e12
        else {
            (p - 1) / (choose(p, num_available_variables) * num_available_variables * (p - num_available_variables))
        }
        })
}


calculate_shap_values <- function(explainer, model, avg_survival_function, data, simplified_inputs, shap_kernel_weights, new_observation, timestamps) {

    w <- shap_kernel_weights
    X <- as.matrix(simplified_inputs)

    R <- solve(t(sqrt(w) * X) %*% (sqrt(w) * X)) %*% t(X * w)

    y <- make_prediction_for_simplified_input(explainer, model, data, simplified_inputs, new_observation, timestamps)

    y <- sweep(y,
               2,
               avg_survival_function)

    R %*% y

}

make_prediction_for_simplified_input <- function(explainer, model, data, simplified_inputs, new_observation, timestamps) {


    preds <- apply(simplified_inputs, 1, function(row) {
        row <- as.logical(row)

        X_tmp <- data
        X_tmp[, row] <- new_observation[, row]
        colnames(X_tmp) <- colnames(data)

        colMeans(explainer$predict_survival_function(model, X_tmp, timestamps))

    })

    return(t(preds))


}

aggregate_surv_shap <- function(survshap, method) {

    switch(method,
           "sum_of_squares" = return(apply(survshap$result, 2, function(x) sum(x^2))),
           "mean_absolute" = return(apply(survshap$result, 2, function(x) mean(abs(x)))),
           "max_absolute" = return(apply(survshap$result, 2, function(x) max(abs(x)))),
           "integral" = return(apply(survshap$result, 2, function(x) {
               x <- abs(x)
               names(x) <- NULL
               times <- survshap$eval_times
               n <- length(x)
               i <- (x[1:(n - 1)] + x[2:n]) * diff(times) / 2
               sum(i) / (max(times) - min(times))
           })),
           stop("aggregation_method has to be one of `sum_of_squares`, `mean_absolute`, `max_absolute` or `integral`"))

}


use_kernelshap <- function(explainer, new_observation, ...){

    predfun <- function(model, newdata){
        explainer$predict_survival_function(
            model,
            newdata,
            times = explainer$times
        )
    }

    tmp_res_list <- sapply(
        X = as.character(seq_len(nrow(new_observation))),
        FUN = function(i) {
            tmp_res <- kernelshap::kernelshap(
                object = explainer$model,
                X = new_observation[as.integer(i), ],
                bg_X = explainer$data,
                pred_fun = predfun,
                verbose = FALSE
            )
            tmp_shap_values <- data.frame(t(sapply(tmp_res$S, cbind)))
            colnames(tmp_shap_values) <- colnames(tmp_res$X)
            rownames(tmp_shap_values) <- paste("t=", explainer$times, sep = "")
            data.table::as.data.table(tmp_shap_values, keep.rownames = TRUE)
        },
        USE.NAMES = TRUE,
        simplify = FALSE
    )

    shap_values <- aggregate_shap_multiple_observations(
        shap_res_list = tmp_res_list,
        feature_names = colnames(new_observation)
    )

    return(shap_values)
}

use_treeshap <- function(explainer, new_observation, ...){

    # init unify_append_args
    unify_append_args <- list()

    if (inherits(explainer$model, "ranger")) {
        # UNIFY_FUN to prepare code for easy Integration of other ml algorithms
        # that are supported by treeshap
        UNIFY_FUN <- treeshap::ranger_surv_fun.unify
        unify_append_args <- list(type = "survival")
    } else {
        stop("Support for `treeshap` is currently only implemented for `ranger`.")
    }

    unify_args <- list(
        rf_model = explainer$model,
        data = explainer$data
    )

    if (length(unify_append_args) > 0) {
        unify_args <- c(unify_args, unify_append_args)
    }

    tmp_unified <- do.call(UNIFY_FUN, unify_args)

    tmp_res_list <- sapply(
        X = as.character(seq_len(nrow(new_observation))),
        FUN = function(i) {
            tmp_res <- do.call(
                rbind,
                lapply(
                    tmp_unified,
                    function(m) {
                    treeshap::treeshap(
                        unified_model = m,
                        x = new_observation
                    )$shaps
                    }
                )
            )

            tmp_shap_values <- data.frame(tmp_res)
            colnames(tmp_shap_values) <- colnames(tmp_res)
            rownames(tmp_shap_values) <- paste("t=", explainer$times, sep = "")
            data.table::as.data.table(tmp_shap_values, keep.rownames = TRUE)
        },
        USE.NAMES = TRUE,
        simplify = FALSE
    )

    shap_values <- aggregate_shap_multiple_observations(
        shap_res_list = tmp_res_list,
        feature_names = colnames(new_observation)
    )

    return(shap_values)
}


aggregate_shap_multiple_observations <- function(shap_res_list, feature_names) {

    if (length(shap_res_list) > 1) {

        full_survshap_results <- data.table::rbindlist(
            l = shap_res_list,
            use.names = TRUE,
            idcol = TRUE
        )

        # compute arithmetic mean for each time-point and feature across
        # multiple observations
        tmp_res <- full_survshap_results[
            , lapply(.SD, mean), by = "rn", .SDcols = feature_names
        ]
    } else {
        # no aggregation required
        tmp_res <- shap_res_list[[1]]
    }
    shap_values <- tmp_res[, .SD, .SDcols = setdiff(colnames(tmp_res), "rn")]
    # transform to data.frame to make everything compatible with
    # previous code
    shap_values <- data.frame(shap_values)
    rownames(shap_values) <- tmp_res$rn
    return(shap_values)
}
