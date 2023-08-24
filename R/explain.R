#' A model-agnostic explainer for survival models
#'
#' Black-box models have vastly different structures. `explain_survival()`
#' returns an explainer object that can be further processed for creating
#' prediction explanations and their visualizations. This function is used to manually
#' create explainers for models not covered by the `survex` package. For selected
#' models the extraction of information can be done automatically. To do
#' this, you can call the `explain()` function for survival models from  `mlr3proba`, `censored`,
#' `randomForestSRC`, `ranger`, `survival` packages and any other model
#' with `pec::predictSurvProb()` method.
#'
#' @param model object - a survival model to be explained
#' @param data data.frame - data which will be used to calculate the explanations. If not provided, then it will be extracted from the model if possible. It should not contain the target columns. NOTE: If the target variable is present in the `data` some functionality breaks.
#' @param y `survival::Surv` object containing event/censoring times and statuses corresponding to `data`
#' @param predict_function  function taking 2 arguments - `model` and `newdata` and returning a single number for each observation - risk score. Observations with higher score are more likely to observe the event sooner.
#' @param predict_function_target_column unused, left for compatibility with DALEX
#' @param residual_function unused, left for compatibility with DALEX
#' @param weights unused, left for compatibility with DALEX
#' @param ... additional arguments, passed to `DALEX::explain()`
#' @param label character - the name of the model. Used to differentiate on visualizations with multiple explainers. By default it's extracted from the 'class' attribute of the model if possible.
#' @param verbose logical, if TRUE (default) then diagnostic messages will be printed
#' @param colorize logical, if TRUE (default) then WARNINGS, ERRORS and NOTES are colorized. Will work only in the R console. By default it is FALSE while knitting and TRUE otherwise.
#' @param model_info a named list (`package`, `version`, `type`) containing information about model. If `NULL`, `survex` will seek for information on its own.
#' @param type type of a model, by default `"survival"`
#'
#' @param times numeric, a vector of times at which the survival function and cumulative hazard function should be evaluated for calculations
#' @param times_generation either `"uniform"` or `"quantiles"`. Sets the way of generating the vector of times based on times provided in the `y` parameter. If `"uniform"` the vector contains 50 equally spaced points between the minimum and maximum observed times; if `"quantiles"` the vector contains 50 points between 0th and 98th percentiles of observed times. Ignored if `times` is not `NULL`.
#' @param predict_survival_function function taking 3 arguments `model`, `newdata` and `times`, and returning a matrix whose each row is a survival function evaluated at `times` for one observation from `newdata`
#' @param predict_cumulative_hazard_function function taking 3 arguments `model`, `newdata` and `times`, and returning a matrix whose each row is a cumulative hazard function evaluated at `times` for one observation from `newdata`
#'
#' @return It is a list containing the following elements:
#'
#' * `model` - the explained model.
#' * `data` - the dataset used for training.
#' * `y` - response for observations from `data`.
#' * `residuals` - calculated residuals.
#' * `predict_function` - function that may be used for model predictions, shall return a single numerical value for each observation.
#' * `residual_function` - function that returns residuals, shall return a single numerical value for each observation.
#' * `class` - class/classes of a model.
#' * `label` - label of explainer.
#' * `model_info` - named list containing basic information about model, like package, version of package and type.
#' * `times` - a vector of times, that are used for evaluation of survival function and cumulative hazard function by default
#' * `predict_survival_function` - function that is used for model predictions in the form of survival function
#' * `predict_cumulative_hazard_function` - function that is used for model predictions in the form of cumulative hazard function
#'
#' @rdname explain_survival
#'
#' @examples
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' cph <- survival::coxph(survival::Surv(time, status) ~ .,
#'     data = veteran,
#'     model = TRUE, x = TRUE
#' )
#' cph_exp <- explain(cph)
#'
#' rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ .,
#'     data = veteran,
#'     respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5
#' )
#' rsf_ranger_exp <- explain(rsf_ranger,
#'     data = veteran[, -c(3, 4)],
#'     y = Surv(veteran$time, veteran$status)
#' )
#'
#' rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' rsf_src_exp <- explain(rsf_src)
#'
#' library(censored, quietly = TRUE)
#'
#' bt <- parsnip::boost_tree() %>%
#'     parsnip::set_engine("mboost") %>%
#'     parsnip::set_mode("censored regression") %>%
#'     generics::fit(survival::Surv(time, status) ~ ., data = veteran)
#' bt_exp <- explain(bt, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status))
#'
#' ###### explain_survival() ######
#'
#' cph <- coxph(Surv(time, status) ~ ., data = veteran)
#'
#' veteran_data <- veteran[, -c(3, 4)]
#' veteran_y <- Surv(veteran$time, veteran$status)
#' risk_pred <- function(model, newdata) predict(model, newdata, type = "risk")
#' surv_pred <- function(model, newdata, times) pec::predictSurvProb(model, newdata, times)
#' chf_pred <- function(model, newdata, times) -log(surv_pred(model, newdata, times))
#'
#' manual_cph_explainer <- explain_survival(
#'     model = cph,
#'     data = veteran_data,
#'     y = veteran_y,
#'     predict_function = risk_pred,
#'     predict_survival_function = surv_pred,
#'     predict_cumulative_hazard_function = chf_pred,
#'     label = "manual coxph"
#' )
#' }
#'
#' @import survival
#' @import ggplot2
#' @import patchwork
#' @importFrom DALEX theme_drwhy theme_drwhy_vertical
#' @importFrom utils tail stack head
#' @importFrom stats median model.frame predict stepfun reorder na.omit aggregate
#'
#' @export
explain_survival <-
    function(model,
             data = NULL,
             y = NULL,
             predict_function = NULL,
             predict_function_target_column = NULL,
             residual_function = NULL,
             weights = NULL,
             ...,
             label = NULL,
             verbose = TRUE,
             colorize = !isTRUE(getOption("knitr.in.progress")),
             model_info = NULL,
             type = NULL,
             times = NULL,
             times_generation = "quantiles",
             predict_survival_function = NULL,
             predict_cumulative_hazard_function = NULL) {
        if (!colorize) {
            color_codes <- list(
                yellow_start = "", yellow_end = "",
                red_start = "", red_end = "",
                green_start = "", green_end = ""
            )
        }


        if (is.null(predict_survival_function) &&
            !is.null(predict_cumulative_hazard_function)) {
            predict_survival_function <- function(model, newdata, times) cumulative_hazard_to_survival(predict_cumulative_hazard_function(model, newdata, times))
            attr(predict_survival_function, "verbose_info") <- "exp(-predict_cumulative_hazard_function) will be used"
            attr(predict_survival_function, "is.default") <- TRUE
        }

        if (is.null(predict_cumulative_hazard_function) &&
            !is.null(predict_survival_function)) {
            predict_cumulative_hazard_function <-
                function(model, newdata, times) survival_to_cumulative_hazard(predict_survival_function(model, newdata, times))
            attr(predict_cumulative_hazard_function, "verbose_info") <- "-log(predict_survival_function) will be used"
            attr(predict_cumulative_hazard_function, "is.default") <- TRUE
        }

        # verbose start
        verbose_cat("Preparation of a new explainer is initiated", verbose = verbose)

        # verbose label
        if (is.null(label)) {
            label <- tail(class(model), 1)
            verbose_cat("  -> model label       : ", label, is.default = TRUE, verbose = verbose)
        } else {
            if (!is.character(label)) {
                label <- substr(as.character(label), 1, 15)
                verbose_cat("  -> model label       : 'label' was not a string class object. Converted. (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("'label' was not a string class object")
            } else if (!is.null(attr(label, "verbose_info")) && attr(label, "verbose_info") == "default") {
                verbose_cat("  -> model label       : ", label, is.default = TRUE, verbose = verbose)
                attr(label, "verbose_info") <- NULL
            } else {
                verbose_cat("  -> model label       : ", label, verbose = verbose)
            }
        }

        # verbose data
        if (is.null(data)) {
            possible_data <- try(model.frame(model), silent = TRUE)
            if (class(possible_data)[1] != "try-error") {
                data <- possible_data
                data <- possible_data[, -1]
                if (is.null(y)) {
                    y <- possible_data[, 1]
                    attr(y, "verbose_info") <- "extracted"
                }
                n <- nrow(data)
                verbose_cat("  -> data              : ", n, " rows ", ncol(data), " cols", "(", color_codes$yellow_start, "extracted from the model", color_codes$yellow_end, ")", verbose = verbose)
            } else {
                # Setting 0 as value of n if data is not present is necessary for future checks
                n <- 0
                verbose_cat("  -> no data available! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("No data available")
            }
        } else {
            n <- nrow(data)
            if (!is.null(attr(data, "verbose_info")) && attr(data, "verbose_info") == "extracted") {
                verbose_cat("  -> data              : ", n, " rows ", ncol(data), " cols", "(", color_codes$yellow_start, "extracted from the model", color_codes$yellow_end, ")", verbose = verbose)
                attr(data, "verbose_info") <- NULL
            } else if (!is.null(attr(data, "verbose_info")) && attr(data, "verbose_info") == "colnames_changed") {
                verbose_cat("  -> data              : ", n, " rows ", ncol(data), " cols", "(", color_codes$yellow_start, "colnames changed to comply with the model", color_codes$yellow_end, ")", verbose = verbose)
                attr(data, "verbose_info") <- NULL
            }
            else {
                verbose_cat("  -> data              : ", n, " rows ", ncol(data), " cols", verbose = verbose)
            }
        }
        if ("tbl" %in% class(data)) {
            data <- as.data.frame(data)
            verbose_cat("  -> data              :  tibble converted into a data.frame", verbose = verbose)
        }

        # verbose target variable
        if (is.null(y)) {
            verbose_cat("  -> target variable   :  not specified! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
            warning("Target variable not specified")
        } else {
            n_events <- sum(y[, 2])
            n_censored <- length(y) - n_events
            frac_censored <- round(n_censored / n, 3)
            if (!is.null(attr(y, "verbose_info")) && attr(y, "verbose_info") == "extracted") {
                verbose_cat("  -> target variable   : ", length(y), " values (", n_events, "events and", n_censored, "censored , censoring rate =", frac_censored, ")", "(", color_codes$yellow_start, "extracted from the model", color_codes$yellow_end, ")", verbose = verbose)
                attr(y, "verbose_info") <- NULL
            } else {
                verbose_cat("  -> target variable   : ", length(y), " values (", n_events, "events and", n_censored, "censored )", verbose = verbose)
            }
            if (length(y) != n) {
                verbose_cat("  -> target variable   :  length of 'y' is different than number of rows in 'data' (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("Length of 'y' is different than number of rows in 'data'")
            }
            if (is.null(data)) {
                verbose_cat("  -> target variable   :  'y' present while 'data' is NULL. (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("'y' present while 'data' is NULL")
            }
        }

        # verbose times
        if (is.null(times)) {
            if (!is.null(y)) {
                switch(times_generation,
                    "uniform" = {
                        times <- seq(min(y[, 1]), max(y[, 1]), length.out = 50)
                        method_description <- "50 uniformly distributed time points from min to max"
                    },
                    "quantiles" = {
                        times <- quantile(y[, 1], seq(0, 0.99, 0.02))
                        method_description <- "50 time points being consecutive quantiles (0.00, 0.02, ..., 0.98)"
                    },
                    stop("times_generation needs to be 'uniform' or 'quantiles'")
                )
                times <- sort(unique(times))
                times_stats <- get_times_stats(times)
                verbose_cat("  -> times             : ", times_stats[1], "unique time points", ", min =", times_stats[2], ", mean =", times_stats[3], ", median =", times_stats[4], ", max =", times_stats[5], verbose = verbose)
                verbose_cat("  -> times             : ", "(", color_codes$yellow_start, paste("generated from y as", method_description), color_codes$yellow_end, ")", verbose = verbose)
            } else {
                verbose_cat("  -> times   :  not specified and automatic generation is impossible ('y' is NULL)! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("'times' not specified and automatic generation is impossible ('y' is NULL)")
            }
        } else {
            times <- sort(unique(times))
            times_stats <- get_times_stats(times)
            verbose_cat("  -> times             : ", times_stats[1], "unique time points", ", min =", times_stats[2], ", mean =", times_stats[3], ", median =", times_stats[4], ", max =", times_stats[5], verbose = verbose)
        }

        # verbose predict function
        if (is.null(predict_function)) {
            if (!is.null(predict_cumulative_hazard_function)) {
                predict_function <- risk_from_chf(predict_cumulative_hazard_function, times)
                verbose_cat("  -> predict function  : ", "sum over the predict_cumulative_hazard_function will be used", is.default = TRUE, verbose = verbose)
            } else {
                verbose_cat("  -> predict function   :  not specified! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("Prediction function not specified")
            }
        } else {
            if (!is.null(attr(predict_function, "verbose_info"))) {
                verbose_cat("  -> predict function  : ", attr(predict_function, "verbose_info"), is.default = attr(predict_function, "is.default"), verbose = verbose)
                attr(predict_function, "verbose_info") <- NULL
                attr(predict_function, "is.default") <- NULL
            } else {
                verbose_cat("  -> predict function  : ", deparse(substitute(predict_function)), verbose = verbose)
            }
            if (!is.null(attr(predict_function, "use.times")) && attr(predict_function, "use.times") == TRUE) {
                predict_function_old <- predict_function
                predict_function <- function(model, newdata) predict_function_old(model, newdata, times = times)
            }
            if (!"function" %in% class(predict_function)) {
                verbose_cat("  -> predict function  :  'predict_function' is not a 'function' class object! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("Prediction function not available")
            }
        }

        # verbose predict survival function
        if (is.null(predict_survival_function)) {
            verbose_cat("  -> predict survival function   :  not specified! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
            warning("Survival function not available")
        } else {
            if (!is.null(attr(predict_survival_function, "verbose_info"))) {
                verbose_cat("  -> predict survival function  : ", attr(predict_survival_function, "verbose_info"), is.default = attr(predict_survival_function, "is.default"), verbose = verbose)
                attr(predict_survival_function, "verbose_info") <- NULL
                attr(predict_survival_function, "is.default") <- NULL
            } else {
                verbose_cat("  -> predict survival function  : ", deparse(substitute(predict_survival_function)), verbose = verbose)
            }
            if (!"function" %in% class(predict_survival_function)) {
                verbose_cat("  -> predict survival function  :  'predict_survival_function' is not a 'function' class object! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("Survival function not available")
            }
        }

        # verbose predict cumulative hazard function
        if (is.null(predict_cumulative_hazard_function)) {
            verbose_cat("  -> predict cumulative hazard function   :  not specified! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
            warning("Cumulative hazard function not available")
        } else {
            if (!is.null(attr(predict_cumulative_hazard_function, "verbose_info"))) {
                verbose_cat("  -> predict cumulative hazard function  : ", attr(predict_cumulative_hazard_function, "verbose_info"), is.default = attr(predict_cumulative_hazard_function, "is.default"), verbose = verbose)
                attr(predict_cumulative_hazard_function, "verbose_info") <- NULL
                attr(predict_cumulative_hazard_function, "is.default") <- NULL
            } else {
                verbose_cat("  -> predict cumulative hazard function  : ", deparse(substitute(predict_cumulative_hazard_function)), verbose = verbose)
            }
            if (!"function" %in% class(predict_cumulative_hazard_function)) {
                verbose_cat("  -> predict cumulative hazard function  :  'predict_cumulative_hazard_function' is not a 'function' class object! (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
                warning("'predict_cumulative_hazard_function' is not a 'function' class object")
            }
        }

        # verbose model info
        if (is.null(model_info)) {
            model_info <- surv_model_info(model)
            verbose_cat("  -> model_info        :  package", model_info$package[1], ", ver.", model_info$ver[1], ", task", model_info$type, is.default = TRUE, verbose = verbose)
        } else {
            verbose_cat("  -> model_info        :  package", model_info$package[1], ", ver.", model_info$ver[1], ", task", model_info$type, verbose = verbose)
        }
        # if type specified then it overwrite the type in model_info
        if (!is.null(type)) {
            model_info$type <- type
            verbose_cat("  -> model_info        :  type set to ", type, verbose = verbose)
        }
        if (class(y)[1] != "Surv") {
            verbose_cat("  -> model_info        :  survival task detected but 'y' is a", class(y)[1], "  (", color_codes$red_start, "WARNING", color_codes$red_end, ")", verbose = verbose)
            verbose_cat("  -> model_info        :  by deafult survival tasks supports only 'y' parameter of 'survival::Surv' class", verbose = verbose)
        }

        explainer <- DALEX::explain(
            model = model,
            data = data,
            y = y,
            predict_function = predict_function,
            predict_function_target_column = NULL,
            residual_function = NULL,
            weights = NULL,
            label = label,
            verbose = FALSE,
            precalculate = FALSE,
            colorize = colorize,
            model_info = model_info,
            type = type,
            times = times,
            predict_survival_function = predict_survival_function,
            predict_cumulative_hazard_function = predict_cumulative_hazard_function,
            ... = ...
        )

        class(explainer) <- c("surv_explainer", class(explainer))

        # verbose end - everything went OK
        verbose_cat("", color_codes$green_start, "A new explainer has been created!", color_codes$green_end, verbose = verbose)
        explainer
    }

#' @rdname explain_survival
#' @export
explain <- function(model,
                    data = NULL,
                    y = NULL,
                    predict_function = NULL,
                    predict_function_target_column = NULL,
                    residual_function = NULL,
                    weights = NULL,
                    ...,
                    label = NULL,
                    verbose = TRUE,
                    colorize = !isTRUE(getOption("knitr.in.progress")),
                    model_info = NULL,
                    type = NULL) {
    UseMethod("explain", model)
}

#' @rdname explain_survival
#' @export
explain.default <- function(model,
                            data = NULL,
                            y = NULL,
                            predict_function = NULL,
                            predict_function_target_column = NULL,
                            residual_function = NULL,
                            weights = NULL,
                            ...,
                            label = NULL,
                            verbose = TRUE,
                            colorize = !isTRUE(getOption("knitr.in.progress")),
                            model_info = NULL,
                            type = NULL) {
    supported_models <- c("aalen", "riskRegression", "cox.aalen", "cph", "coxph", "selectCox", "pecCforest", "prodlim", "psm", "survfit", "pecRpart")
    if (inherits(model, supported_models)) {
        return(
            explain_survival(
                model,
                data = data,
                y = y,
                predict_function = predict_function,
                predict_function_target_column = predict_function_target_column,
                residual_function = residual_function,
                weights = weights,
                ...,
                label = label,
                verbose = verbose,
                colorize = colorize,
                model_info = model_info,
                type = type,
                predict_survival_function = pec::predictSurvProb
            )
        )
    }
    if (inherits(model, "sksurv.base.SurvivalAnalysisMixin")){
        return(
            explain.sksurv(model,
                data = data,
                y = y,
                predict_function = predict_function,
                predict_function_target_column = predict_function_target_column,
                residual_function = residual_function,
                weights = weights,
                ...,
                label = label,
                verbose = verbose,
                colorize = colorize,
                model_info = model_info,
                type = type
            )
        )
    }

    DALEX::explain(model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = !isTRUE(getOption("knitr.in.progress")),
        model_info = model_info,
        type = type
    )

}

#' @export
explain.coxph <- function(model,
                          data = NULL,
                          y = NULL,
                          predict_function = NULL,
                          predict_function_target_column = NULL,
                          residual_function = NULL,
                          weights = NULL,
                          ...,
                          label = NULL,
                          verbose = TRUE,
                          colorize = !isTRUE(getOption("knitr.in.progress")),
                          model_info = NULL,
                          type = NULL,
                          times = NULL,
                          times_generation = "quantiles",
                          predict_survival_function = NULL,
                          predict_cumulative_hazard_function = NULL) {
    if (is.null(data)) {
        data <- model$model[, attr(model$terms, "term.labels")]
        if (is.null(data)) {
            stop(
                "use `model=TRUE` and `x=TRUE` while creating coxph model or provide `data` manually"
            )
        }
        attr(data, "verbose_info") <- "extracted"
    }

    if (is.null(y)) {
        y <- model$y
        if (is.null(y)) {
            stop("use `y=TRUE` while creating coxph model or provide `y` manually")
        }
        attr(y, "verbose_info") <- "extracted"
    }

    if (is.null(predict_survival_function)) {
        predict_survival_function <- function(model, newdata, times) {
            pec::predictSurvProb(model, newdata, times)
        }
        attr(predict_survival_function, "verbose_info") <- "predictSurvProb.coxph will be used"
        attr(predict_survival_function, "is.default") <- TRUE
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        predict_cumulative_hazard_function <-
            function(model, newdata, times) {
                survival_to_cumulative_hazard(predict_survival_function(model, newdata, times))
            }
        attr(predict_cumulative_hazard_function, "verbose_info") <- "-log(predict_survival_function) will be used"
        attr(predict_cumulative_hazard_function, "is.default") <- TRUE
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        predict_function <- function(model, newdata) {
            predict(model, newdata, type = "risk")
        }
        attr(predict_function, "verbose_info") <- "predict.coxph with type = 'risk' will be used"
        attr(predict_function, "is.default") <- TRUE
    } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}


#' @export
explain.ranger <- function(model,
                           data = NULL,
                           y = NULL,
                           predict_function = NULL,
                           predict_function_target_column = NULL,
                           residual_function = NULL,
                           weights = NULL,
                           ...,
                           label = NULL,
                           verbose = TRUE,
                           colorize = !isTRUE(getOption("knitr.in.progress")),
                           model_info = NULL,
                           type = NULL,
                           times = NULL,
                           times_generation = "quantiles",
                           predict_survival_function = NULL,
                           predict_cumulative_hazard_function = NULL) {
    if (is.null(predict_survival_function)) {
        predict_survival_function <- transform_to_stepfunction(predict,
            type = "survival",
            times_element = "unique.death.times",
            prediction_element = "survival"
        )
        attr(predict_survival_function, "verbose_info") <- "stepfun based on predict.ranger()$survival will be used"
        attr(predict_survival_function, "is.default") <- TRUE
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        predict_cumulative_hazard_function <- transform_to_stepfunction(predict,
            type = "chf",
            times_element = "unique.death.times",
            prediction_element = "chf"
        )
        attr(predict_cumulative_hazard_function, "verbose_info") <- "stepfun based on predict.ranger()$chf will be used"
        attr(predict_cumulative_hazard_function, "is.default") <- TRUE
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        predict_function <- function(model, newdata, times) {
            rowSums(predict_cumulative_hazard_function(model, newdata, times))
        }
        attr(predict_function, "verbose_info") <- "sum over the predict_cumulative_hazard_function will be used"
        attr(predict_function, "is.default") <- TRUE
        attr(predict_function, "use.times") <- TRUE
    } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}


#' @export
explain.rfsrc <- function(model,
                          data = NULL,
                          y = NULL,
                          predict_function = NULL,
                          predict_function_target_column = NULL,
                          residual_function = NULL,
                          weights = NULL,
                          ...,
                          label = NULL,
                          verbose = TRUE,
                          colorize = !isTRUE(getOption("knitr.in.progress")),
                          model_info = NULL,
                          type = NULL,
                          times = NULL,
                          times_generation = "quantiles",
                          predict_survival_function = NULL,
                          predict_cumulative_hazard_function = NULL) {
    if (is.null(label)) {
        label <- class(model)[1]
        attr(label, "verbose_info") <- "default"
    }

    if (is.null(data)) {
        data <- model$xvar
        attr(data, "verbose_info") <- "extracted"
    }

    if (is.null(y)) {
        tmp_y <- model$yvar
        y <- survival::Surv(tmp_y[, 1], tmp_y[, 2])
        attr(y, "verbose_info") <- "extracted"
    }

    if (is.null(predict_survival_function)) {
        predict_survival_function <- transform_to_stepfunction(predict,
            type = "survival",
            times_element = "time.interest",
            prediction_element = "survival"
        )
        attr(predict_survival_function, "verbose_info") <- "stepfun based on predict.rfsrc()$survival will be used"
        attr(predict_survival_function, "is.default") <- TRUE
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        predict_cumulative_hazard_function <- transform_to_stepfunction(predict,
            type = "chf",
            times_element = "time.interest",
            prediction_element = "chf"
        )
        attr(predict_cumulative_hazard_function, "verbose_info") <- "stepfun based on predict.rfsrc()$chf will be used"
        attr(predict_cumulative_hazard_function, "is.default") <- TRUE
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        predict_function <- function(model, newdata, times) {
            rowSums(predict_cumulative_hazard_function(model, newdata, times = times))
        }
        attr(predict_function, "verbose_info") <- "sum over the predict_cumulative_hazard_function will be used"
        attr(predict_function, "is.default") <- TRUE
        attr(predict_function, "use.times") <- TRUE
    } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}


#' @export
explain.model_fit <- function(model,
                              data = NULL,
                              y = NULL,
                              predict_function = NULL,
                              predict_function_target_column = NULL,
                              residual_function = NULL,
                              weights = NULL,
                              ...,
                              label = NULL,
                              verbose = TRUE,
                              colorize = !isTRUE(getOption("knitr.in.progress")),
                              model_info = NULL,
                              type = NULL,
                              times = NULL,
                              times_generation = "quantiles",
                              predict_survival_function = NULL,
                              predict_cumulative_hazard_function = NULL) {
    if (is.null(label)) {
        label <- paste(rev(class(model)), collapse = "")
        attr(label, "verbose_info") <- "default"
    }

    if (is.null(predict_survival_function)) {
        predict_survival_function <- function(model, newdata, times) {
            prediction <- predict(model, new_data = newdata, type = "survival", eval_time = times)$.pred
            return_matrix <- t(sapply(prediction, function(x) x$.pred_survival))
            return_matrix[is.na(return_matrix)] <- 0
            return_matrix
        }
        attr(predict_survival_function, "verbose_info") <- "predict.model_fit with type = 'survival' will be used"
        attr(predict_survival_function, "is.default") <- TRUE
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        predict_cumulative_hazard_function <-
            function(object, newdata, times) survival_to_cumulative_hazard(predict_survival_function(object, newdata, times))
        attr(predict_cumulative_hazard_function, "verbose_info") <- "-log(predict_survival_function) will be used"
        attr(predict_cumulative_hazard_function, "is.default") <- TRUE
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        if (model$spec$engine %in% c("mboost", "survival", "glmnet", "flexsurv", "flexsurvspline")) {
            predict_function <- function(model, newdata, times) predict(model, new_data = newdata, type = "linear_pred")$.pred_linear_pred
            attr(predict_function, "verbose_info") <- "predict.model_fit with type = 'linear_pred' will be used"
        } else {
            predict_function <- function(model, newdata, times) rowSums(predict_cumulative_hazard_function(model, newdata, times = times))
            attr(predict_function, "verbose_info") <- "sum over the predict_cumulative_hazard_function will be used"
        }
        attr(predict_function, "use.times") <- TRUE
        attr(predict_function, "is.default") <- TRUE
    } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}


#' @export
explain.LearnerSurv <- function(model,
                                data = NULL,
                                y = NULL,
                                predict_function = NULL,
                                predict_function_target_column = NULL,
                                residual_function = NULL,
                                weights = NULL,
                                ...,
                                label = NULL,
                                verbose = TRUE,
                                colorize = !isTRUE(getOption("knitr.in.progress")),
                                model_info = NULL,
                                type = NULL,
                                times = NULL,
                                times_generation = "quantiles",
                                predict_survival_function = NULL,
                                predict_cumulative_hazard_function = NULL) {
    if (is.null(label)) {
        label <- class(model)[1]
        attr(label, "verbose_info") <- "default"
    }

    if (is.null(predict_survival_function)) {
        if ("distr" %in% model$predict_types) {
            predict_survival_function <- function(model, newdata, times) t(model$predict_newdata(newdata)$distr$survival(times))
            attr(predict_survival_function, "verbose_info") <- "predict_newdata()$distr$survival will be used"
            attr(predict_survival_function, "is.default") <- TRUE
        }
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        if ("distr" %in% model$predict_types) {
            predict_cumulative_hazard_function <- function(model, newdata, times) t(model$predict_newdata(newdata)$distr$cumHazard(times))
            attr(predict_cumulative_hazard_function, "verbose_info") <- "predict_newdata()$distr$cumHazard will be used"
            attr(predict_cumulative_hazard_function, "is.default") <- TRUE
        }
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        if ("crank" %in% model$predict_types) {
            predict_function <- function(model, newdata, times) model$predict_newdata(newdata)$crank
            attr(predict_function, "verbose_info") <- "predict_newdata()$crank will be used"
        } else {
            predict_function <- function(model, newdata, times) {
                rowSums(predict_cumulative_hazard_function(model, newdata, times))
            }
            attr(predict_function, "verbose_info") <- "sum over the predict_cumulative_hazard_function will be used"
        }
        attr(predict_function, "is.default") <- TRUE
        attr(predict_function, "use.times") <- TRUE
    } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}

#' @export
explain.sksurv <- function(model,
                          data = NULL,
                          y = NULL,
                          predict_function = NULL,
                          predict_function_target_column = NULL,
                          residual_function = NULL,
                          weights = NULL,
                          ...,
                          label = NULL,
                          verbose = TRUE,
                          colorize = !isTRUE(getOption("knitr.in.progress")),
                          model_info = NULL,
                          type = NULL,
                          times = NULL,
                          times_generation = "quantiles",
                          predict_survival_function = NULL,
                          predict_cumulative_hazard_function = NULL){
    if (is.null(label)) {
        label <- class(model)[1]
        attr(label, "verbose_info") <- "default"
    }

    if (is.null(predict_survival_function)) {
        if (reticulate::py_has_attr(model, "predict_survival_function")) {
            predict_survival_function <- function(model, newdata, times){
                raw_preds <- model$predict_survival_function(newdata)
                t(sapply(raw_preds, function(sf) as.vector(sf(times))))
            }
            attr(predict_survival_function, "verbose_info") <- "predict_survival_function from scikit-survival will be used"
            attr(predict_survival_function, "is.default") <- TRUE
        }
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        if (reticulate::py_has_attr(model, "predict_cumulative_hazard_function")) {
            predict_cumulative_hazard_function <- function(model, newdata, times){
                raw_preds <- model$predict_cumulative_hazard_function(newdata)
                t(sapply(raw_preds, function(chf) as.vector(chf(times))))
            }
            attr(predict_cumulative_hazard_function, "verbose_info") <- "predict_cumulative_hazard_function from scikit-survival will be used"
            attr(predict_cumulative_hazard_function, "is.default") <- TRUE
        }
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        if (reticulate::py_has_attr(model, "predict")) {
            predict_function <- function(model, newdata) model$predict(newdata)
            attr(predict_function, "verbose_info") <- "predict from scikit-survival will be used"
            attr(predict_function, "is.default") <- TRUE
        } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
        }
    }

    if (!is.null(data) & any(colnames(data) != model$feature_names_in_)) {
        colnames(data) <- sub("[.]", "=", colnames(data))
        attr(data, "verbose_info") <- "colnames_changed"
    }

    if (class(model)[1] != "sksurv"){
        class(model) <- c("sksurv", class(model))
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}


#' @export
explain.flexsurvreg <- function(model,
                                data = NULL,
                                y = NULL,
                                predict_function = NULL,
                                predict_function_target_column = NULL,
                                residual_function = NULL,
                                weights = NULL,
                                ...,
                                label = NULL,
                                verbose = TRUE,
                                colorize = !isTRUE(getOption("knitr.in.progress")),
                                model_info = NULL,
                                type = NULL,
                                times = NULL,
                                times_generation = "quantiles",
                                predict_survival_function = NULL,
                                predict_cumulative_hazard_function = NULL) {
    if (is.null(label)) {
        label <- class(model)[1]
        attr(label, "verbose_info") <- "default"
    }

    if (is.null(predict_survival_function)) {
        predict_survival_function <-  function(model, newdata, times){
                    raw_preds <- predict(model, newdata = newdata, times = times, type = "survival")
                    preds <- do.call(rbind, lapply(raw_preds[[1]], function(x) t(x[".pred_survival"])))
                    rownames(preds) <- NULL
                    preds
                }
        attr(predict_survival_function, "verbose_info") <- "predict.flexsurvreg with type = 'survival' will be used"
        attr(predict_survival_function, "is.default") <- TRUE
    } else {
        attr(predict_survival_function, "verbose_info") <- deparse(substitute(predict_survival_function))
    }

    if (is.null(predict_cumulative_hazard_function)) {
        predict_cumulative_hazard_function <- function(model, newdata, times){
            raw_preds <- predict(model, newdata = newdata, times = times, type = "cumhaz")
            preds <- do.call(rbind, lapply(raw_preds[[1]], function(x) t(x[".pred_cumhaz"])))
            rownames(preds) <- NULL
            preds
        }
        attr(predict_cumulative_hazard_function, "verbose_info") <- "predict.flexsurvreg with type = 'cumhaz' will be used"
        attr(predict_cumulative_hazard_function, "is.default") <- TRUE
    } else {
        attr(predict_cumulative_hazard_function, "verbose_info") <- deparse(substitute(predict_cumulative_hazard_function))
    }

    if (is.null(predict_function)) {
        predict_function <- function(model, newdata){
            predict(model, newdata = newdata, type = "link")[[".pred_link"]]
        }
        attr(predict_function, "verbose_info") <- "predict.flexsurvreg with type = 'link' will be used"
        attr(predict_function, "is.default") <- TRUE
    } else {
        attr(predict_function, "verbose_info") <- deparse(substitute(predict_function))
    }

    possible_data <- model.frame(parmodel)
    if (is.null(data)) {
        data <- possible_data[,-c(1, ncol(possible_data))]
        attr(data, "verbose_info") <- "extracted"
    }

    if (is.null(y)) {
        y <- possible_data[,1]
        attr(y, "verbose_info") <- "extracted"
    }

    explain_survival(
        model,
        data = data,
        y = y,
        predict_function = predict_function,
        predict_function_target_column = predict_function_target_column,
        residual_function = residual_function,
        weights = weights,
        ... = ...,
        label = label,
        verbose = verbose,
        colorize = colorize,
        model_info = model_info,
        type = type,
        times = times,
        times_generation = times_generation,
        predict_survival_function = predict_survival_function,
        predict_cumulative_hazard_function = predict_cumulative_hazard_function
    )
}

verbose_cat <- function(..., is.default = NULL, verbose = TRUE) {
    if (verbose) {
        if (!is.null(is.default)) {
            txt <- paste(..., "(", color_codes$yellow_start, "default", color_codes$yellow_end, ")")
            cat(txt, "\n")
        } else {
            cat(..., "\n")
        }
    }
}

get_times_stats <- function(times) {
    c(length(times), min(times), mean(times), median(times), max(times))
}

#
# colors for WARNING, NOTE, DEFAULT
#
color_codes <- list(
    yellow_start = "\033[33m", yellow_end = "\033[39m",
    red_start = "\033[31m", red_end = "\033[39m",
    green_start = "\033[32m", green_end = "\033[39m"
)
