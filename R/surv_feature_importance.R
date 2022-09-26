#' Helper functions for `model_parts.R`
#' @rdname surv_feature_importance
#' @keywords internal
surv_feature_importance <- function(x, ...) UseMethod("surv_feature_importance", x)


#' Helper functions for `model_parts.R`
#'
#' @rdname surv_feature_importance
#'
#' @param x an explainer object - model preprocessed by the `explain()` function
#' @param loss_function a function that will be used to assess variable importance, by default `loss_brier_score` for survival models. The function can be supplied manually but has to have these named parameters (`y_true`, `risk`, `surv`, `times`), where `y_true` represents the `survival::Surv` object with observed times and statuses, `risk` is the risk score calculated by the model, and `surv` is the survival function for each observation evaluated at `times`.
#' @param ... other parameters, currently ignored
#' @param type a character vector, if `"raw"` the results are losses after the permutation, if `"ratio"` the results are in the form `loss/loss_full_model` and if `"difference"` the results are of the form loss - loss_full_model`
#' @param B numeric, number of permutations to be calculated
#' @param variables a character vector, names of variables to be included in the calculation
#' @param variable_groups a list of character vectors of names of explanatory variables. For each vector, a single variable importance measure is computed for the joint effect of the variables which names are provided in the vector. By default, `variable_groups = NULL`, in which case variable importance measures are computed separately for all variables indicated in the `variables` argument.
#' @param N numeric, number of observations that are to be sampled from the dataset for the purpose of calculation
#' @param label label of the model, if provided overrides `explainer$label`
#'
#' @return A data.frame containing results of the calculation.
#'
#' @details
#' *Note*: This function can be run within `progressr::with_progress()` to display a progress bar, as the execution can take long, especially on large datasets.
#'
#' @keywords internal
surv_feature_importance.surv_explainer <- function(x,
                                         loss_function = NULL,
                                         ...,
                                         type = c("raw", "ratio", "difference"),
                                         B = 10,
                                         variables = NULL,
                                         variable_groups = NULL,
                                         N = NULL,
                                         label = NULL) {
    if (is.null(x$data)) stop("The feature_importance() function requires explainers created with specified 'data' parameter.")
    if (is.null(x$y)) stop("The feature_importance() function requires explainers created with specified 'y' parameter.")
    if (is.null(x$predict_survival_function)) stop("The feature_importance() function requires explainers created with specified 'predict_survival_function' parameter.")



    model <- x$model
    data <- x$data
    predict_function <- x$predict_function
    predict_survival_function <- x$predict_survival_function
    if (is.null(label)) {
        label <- x$label
    }
    y <- x$y
    times <- x$times


    surv_feature_importance.default(model,
                               data,
                               y,
                               times,
                               predict_function = predict_function,
                               predict_survival_function = predict_survival_function,
                               loss_function = loss_function,
                               label = label,
                               type = type,
                               N = N,
                               B = B,
                               variables = variables,
                               variable_groups = variable_groups,
                               ...
    )

}


surv_feature_importance.default <- function(x,
                                                       data,
                                                       y,
                                                       times,
                                                       predict_function = NULL,
                                                       predict_survival_function = NULL,
                                                       loss_function = DALEX::loss_root_mean_square,
                                                       ...,
                                                       label = class(x)[1],
                                                       type = c("raw", "ratio", "difference"),
                                                       B = 10,
                                                       variables = NULL,
                                                       N = NULL,
                                                       variable_groups = NULL) {




    if (!is.null(variable_groups)) {
        if (!inherits(variable_groups, "list")) stop("variable_groups should be of class list")

        wrong_names <- !all(sapply(variable_groups, function(variable_set) {
            all(variable_set %in% colnames(data))
        }))

        if (wrong_names) stop("You have passed wrong variables names in variable_groups argument")
        if (!all(sapply(variable_groups, class) == "character")) stop("Elements of variable_groups argument should be of class character")
        if (is.null(names(variable_groups))) warning("You have passed an unnamed list. The names of variable groupings will be created from variables names.")
    }

    type <- match.arg(type)
    B <- max(1, round(B))

    # Adding names for variable_groups if not specified
    if (!is.null(variable_groups) && is.null(names(variable_groups))) {
        names(variable_groups) <- sapply(variable_groups, function(variable_set) {
            paste0(variable_set, collapse = "; ")
        })
    }

    # if `variable_groups` are not specified, then extract from `variables`
    if (is.null(variable_groups)) {
        # if `variables` are not specified, then extract from data
        if (is.null(variables)) {
            variables <- colnames(data)
            names(variables) <- colnames(data)
        }
    } else {
        variables <- variable_groups
    }

    # start: actual calculations
    # one permutation round: subsample data, permute variables and compute losses
    if (requireNamespace("progressr", quietly = TRUE)) prog <- progressr::progressor(along = 1:((length(variables) + 2) * B))
    sampled_rows <- 1:nrow(data)
    loss_after_permutation <- function() {
        if (!is.null(N)) {
            if (N < nrow(data)) {
                # sample N points
                sampled_rows <- sample(1:nrow(data), N)
            }
        }
        sampled_data <- data[sampled_rows, , drop = FALSE]
        observed <- y[sampled_rows]

        surv_true <- predict_survival_function(x, sampled_data, times)

        risk_true <- predict_function(x, sampled_data)
        # loss on the full model or when outcomes are permuted
        loss_full <- loss_function(observed, risk_true, surv_true, times)
        if (requireNamespace("progressr", quietly = TRUE)) prog()
        chosen <- sample(1:nrow(observed))
        loss_baseline <- loss_function(observed[chosen, ],  risk_true, surv_true, times)
        if (requireNamespace("progressr", quietly = TRUE)) prog()
        # loss upon dropping a single variable (or a single group)
        loss_variables <- sapply(variables, function(variables_set) {
            ndf <- sampled_data
            ndf[, variables_set] <- ndf[sample(1:nrow(ndf)), variables_set]
            predicted <- predict_function(x, ndf)
            predicted_surv <- predict_survival_function(x, ndf, times)
            if (requireNamespace("progressr", quietly = TRUE)) prog()
            loss_function(observed, predicted, predicted_surv, times)
        })

        cbind(data.frame(times = times), "_full_model_" = loss_full, loss_variables, "_baseline_" = loss_baseline)
    }
    # permute B times, collect results into single matrix
    raw <- do.call("rbind", replicate(B, loss_after_permutation(), simplify = FALSE))
    raw$permutation <- rep(1:B, each = length(times))
    tmp <- aggregate(. ~ times, raw, mean)
    tmp$permutation <- rep(0, times = nrow(tmp))

    res <- rbind(tmp, raw)
    res$label <- rep(label, times = nrow(res))

    if (type %in% c("ratio", "difference")) {
        res_full <- res[res$permutation == 0, c("times", "_full_model_")]
        colnames(res_full) <- c("times", "reference")
        res <- merge(res, res_full, by = "times")

    }
    if (type == "ratio") {

        res[, 2:(ncol(res) - 3)] <- res[, 2:(ncol(res) - 3)] / res[, (ncol(res))]
        res$reference <- NULL

    }
    if (type == "difference") {
        res[, 2:(ncol(res) - 3)] <- res[, 2:(ncol(res) - 3)] - res[, (ncol(res))]
        res$reference <- NULL
    }

    # record details of permutations
    attr(res, "B") <- B

    ret <- list(result = res, eval_times = unique(res$times))
    class(ret) <- c("surv_feature_importance", "list")

    if (!is.null(attr(loss_function, "loss_name"))) {
        attr(ret, "loss_name") <- attr(loss_function, "loss_name")
    }

    attr(ret, "type") <- type

    ret
}
