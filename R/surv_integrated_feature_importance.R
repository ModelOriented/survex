#' Helper functions for `model_parts.R`
#'
#' This function is used to calculate permutational feature importance using
#' an aggregate (for now only integral) of a time dependent metric. The result is
#' the combined change in loss function across all time points - a single value.
#'
#'
#' @param x an explainer object - model preprocessed by the `explain()` function
#' @param loss_function a function that will be used to assess variable importance, by default `loss_brier_score` for survival models. The function should take can be supplied manually but has to have these named parameters (`y_true`, `risk`, `surv`, `times`), where `y_true` represents the `survival::Surv` object with observed times and statuses, `risk` is the risk score calculated by the model, and `surv` is the survival function for each observation evaluated at `times`
#' @param ... other parameters, currently ignored
#' @param type a character vector, if `"raw"` the results are losses after the permutation, if `"ratio"` the results are in the form `loss/loss_full_model` and if `"difference"` the results are of the form loss - loss_full_model`
#' @param B numeric, number of permutations to be calculated
#' @param variables a character vector, names of variables to be included in the calculation
#' @param variable_groups a list of character vectors of names of explanatory variables. For each vector, a single variable-importance measure is computed for the joint effect of the variables which names are provided in the vector. By default, variable_groups = NULL, in which case variable-importance measures are computed separately for all variables indicated in the variables argument
#' @param N numeric, number of observations that are to be sampled from the dataset for the purpose of calculation
#' @param label label of the model, if provides overrides x$label
#'
#' @return A data.frame containing results of the calculation.
#'
#' @details
#' *Note*: This function can be run within `progressr::with_progress()` to display a progress bar, as the execution can take long, especially on large datasets.
#'
#' @keywords internal
surv_integrated_feature_importance <- function(x,
                                               loss_function = DALEX::loss_root_mean_square,
                                               ...,
                                               type = c("raw", "ratio", "difference"),
                                               B = 10,
                                               variables = NULL,
                                               variable_groups = NULL,
                                               N = NULL,
                                               label = NULL) {

    if (is.null(x$data)) stop("The feature_importance() function requires explainers created with specified 'data' parameter.")
    if (is.null(x$y)) stop("The feature_importance() function requires explainers created with specified 'y' parameter.")

    # extracts model, data and predict function from the explainer
    explainer <- x
    data <- x$data
    predict_function <- x$predict_function
    predict_survival_function <- x$predict_survival_function
    if (is.null(label)) {
        label <- x$label
    }
    y <- x$y
    x <- x$model
    times <- explainer$times


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
    prog <- progressr::progressor(along = 1:((length(variables) + 2) * B))
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
        # loss on the full model or when outcomes are permuted

        risks <- predict_function(x, sampled_data)
        survs <- predict_survival_function(x, sampled_data, times)


        loss_full <- loss_function(observed, risks, survs, times)
        prog()
        loss_baseline <- loss_function(sample(observed), risks, survs, times)
        prog()
        # loss upon dropping a single variable (or a single group)
        loss_variables <- sapply(variables, function(variables_set) {
            ndf <- sampled_data
            ndf[, variables_set] <- ndf[sample(1:nrow(ndf)), variables_set]
            predicted_risks <- predict_function(x, ndf)
            predicted_survs <- predict_survival_function(x, ndf, times)
            prog()
            loss_function(observed, predicted_risks, predicted_survs, times)

        })
        c("_full_model_" = loss_full, loss_variables, "_baseline_" = loss_baseline)
    }
    # permute B times, collect results into single matrix
    raw <- replicate(B, loss_after_permutation())

    # main result df with dropout_loss averages, with _full_model_ first and _baseline_ last
    res <- apply(raw, 1, mean)
    res_baseline <- res["_baseline_"]
    res_full <- res["_full_model_"]
    res <- sort(res[!names(res) %in% c("_full_model_", "_baseline_")])
    res <- data.frame(
        variable = c("_full_model_", names(res), "_baseline_"),
        permutation = 0,
        dropout_loss = c(res_full, res, res_baseline),
        label = label,
        row.names = NULL
    )
    if (type == "ratio") {
        res$dropout_loss = res$dropout_loss / res_full
    }
    if (type == "difference") {
        res$dropout_loss = res$dropout_loss - res_full
    }


    # record details of permutations
    attr(res, "B") <- B

    if (B > 1) {
        res_B <- data.frame(
            variable = rep(rownames(raw), ncol(raw)),
            permutation = rep(seq_len(B), each = nrow(raw)),
            dropout_loss = as.vector(raw),
            label = label
        )

        # here mean full model is used (full model for given permutation is an option)
        if (type == "ratio") {
            res_B$dropout_loss = res_B$dropout_loss / res_full
        }
        if (type == "difference") {
            res_B$dropout_loss = res_B$dropout_loss - res_full
        }

        res <- rbind(res, res_B)
    }

    class(res) <- c("feature_importance_explainer", "data.frame")

    if (!is.null(attr(loss_function, "loss_name"))) {
        attr(res, "loss_name") <- attr(loss_function, "loss_name")
    }
    res

}
