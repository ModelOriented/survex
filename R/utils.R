#' Transform Cumulative Hazard to Survival
#'
#' Helper function to transform between CHF and survival function
#'
#' @param hazard_functions matrix or vector, with each row representing a cumulative hazard function
#'
#' @return A matrix or vector transformed to the form of a survival function.
#'
#' @examples
#' library(survex)
#'
#' vec <- c(1, 2, 3, 4, 5)
#' matr <- matrix(c(1, 2, 3, 2, 4, 6), ncol = 3)
#'
#' cumulative_hazard_to_survival(vec)
#'
#' cumulative_hazard_to_survival(matr)
#'
#' @export
cumulative_hazard_to_survival <- function(hazard_functions) {
        return(exp(-hazard_functions))
}

#' Transform Survival to Cumulative Hazard
#'
#' Helper function to transform between survival function and CHF
#' @param survival_functions matrix or vector, with each row representing a survival function
#' @param epsilon a positive numeric number to add, so that the logarithm can be taken
#'
#' @return A matrix or vector transformed to the form of a cumulative hazard function.
#'
#' @examples
#' library(survex)
#'
#' vec <- c(1, 0.9, 0.8, 0.7, 0.6)
#' matr <- matrix(c(1, 0.9, 0.8, 1, 0.8, 0.6), ncol = 3)
#'
#' survival_to_cumulative_hazard(vec)
#'
#' survival_to_cumulative_hazard(matr)
#'
#' @export
survival_to_cumulative_hazard <- function(survival_functions, epsilon = 0) {
        return(-log(survival_functions))
}

# tests if the explainer has all the required fields
test_explainer <- function(explainer,
                           function_name,
                           has_data = FALSE,
                           has_y = FALSE,
                           has_survival = FALSE,
                           has_chf = FALSE,
                           has_predict = FALSE)
{
    if (!("surv_explainer" %in% class(explainer)))
       stop(paste0("The ", function_name, " function requires an object created with survex::explain() function."))
    if (has_data && is.null(explainer$data))
       stop(paste0("The ", function_name, " function requires explainers with specified `data` parameter"))
    if (has_y && is.null(explainer$y))
        stop(paste0("The ", function_name, " function requires explainers with specified `y` parameter"))
    if (has_survival && is.null(explainer$predict_survival_function))
        stop(paste0("The ", function_name, " function requires explainers with specified `predict_survival_function` parameter"))
    if (has_chf && is.null(explainer$predict_cumulative_hazard_function))
        stop(paste0("The ", function_name, " function requires explainers with specified `predict_cumulative_hazard_function` parameter"))
    if (has_predict && is.null(explainer$predict_function))
        stop(paste0("The ", function_name, " function requires explainers with specified `predict_risk` parameter"))
}


#' @importFrom DALEX colors_discrete_drwhy
generate_discrete_color_scale <- function(n, colors = NULL) {

    if (is.null(colors) || length(colors) < n) return(colors_discrete_drwhy(n))
    else return(colors[(0:(n - 1) %% length(colors)) + 1])

}

#' Transform Fixed Point Prediction into a Stepfunction
#'
#' Some models return the survival function or cumulative hazard function prediction at the times of events present in the training data set. This is a convenient utility to allow the prediction to be evaluated at any time.
#'
#' @param predict_function a function making the prediction based on `model` and `newdata` arguments, the `...` parameter is also passed to this function. It has to return either a numeric vector of the same length as `eval_times`, a matrix with this number of columns and the same number of rows as `nrow(newdata)`. It can also return a list, with one of the elements containing such an object.
#' @param eval_times a numeric vector of times, at which the fixed predictions are made. This can be `NULL`, if `predict_function` returns a list which contains such a vector.
#' @param ... other parameters passed to predict_function
#' @param type the type of function to be returned, either `"survival"`, `"chf"` or `NULL` this chooses the value of the step function before the first prediction time. If `"survival"` then it is 1, if `"chf"` then 0, otherwise, it is the value of the prediction for the first time in numerical order.
#' @param prediction_element if `predict_function` returns a list with the matrix as one of its elements, this parameter should contain the name of this element
#' @param times_element if `predict_function` returns a list with the matrix as one of its elements, this parameter should contain the name of this element
#'
#' @return The function returns a function with three arguments, (`model`, `newdata`, `times`), ready to supply it to an explainer.
#'
#' @examples
#' library(survex)
#' library(survival)
#'
#' rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#'
#' chf_function <- transform_to_stepfunction(predict,
#'                                           type = "chf",
#'                                           prediction_element = "chf",
#'                                           times_element = "time.interest")
#'
#' explainer <- explain(rsf_src, predict_cumulative_hazard_function = chf_function)
#'
#' @export
transform_to_stepfunction <- function(predict_function, eval_times = NULL, ..., type = NULL, prediction_element = NULL, times_element = NULL) {

        function(model, newdata, times) {
            raw_prediction <- predict_function(model, newdata, ...)
            if (!is.null(times_element)) eval_times <- raw_prediction[[times_element]]
            if (!is.null(prediction_element)) prediction <- raw_prediction[[prediction_element]]
            n_rows <- ifelse(is.null(dim(prediction)), 1, nrow(prediction))
            return_matrix <- matrix(nrow = n_rows, ncol = length(times))


            if (is.null(dim(prediction))) {
                padding <- switch(type,
                                  "survival" = 1,
                                  "chf" = 0,
                                  prediction[1])
                stepfunction <- stepfun(eval_times, c(padding, prediction))
                return_matrix[1, ] <- stepfunction(times)

            }
            else {
                for (i in 1:n_rows) {
                    padding <- switch(type,
                                      "survival" = 1,
                                      "chf" = 0,
                                      prediction[i, 1])
                    stepfunction <- stepfun(eval_times, c(padding, prediction[i, ]))
                    return_matrix[i, ] <- stepfunction(times)
                }
            }

            return_matrix
        }

}

#' Generate Risk Prediction based on the Survival Function
#'
#' Some models do not come with a ready to use risk prediction. This function allows for its generation based on the cumulative hazard function.
#'
#' @param predict_cumulative_hazard_function a function of three arguments (`model`, `newdata`, `times`) that allows for making cumulative hazard predictions.
#' @param times a numeric vector of times at which the function should be evaluated.
#'
#' @return A function of two arguments (`model`, `newdata`) returning a vector of risks.
#'
#' @examples
#' library(survex)
#' library(survival)
#'
#' rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#'
#' chf_function <- transform_to_stepfunction(predict,
#'                                           type = "chf",
#'                                           prediction_element = "chf",
#'                                           times_element = "time.interest")
#' risk_function <- risk_from_chf(chf_function, unique(veteran$time))
#'
#' explainer <- explain(rsf_src,
#'                      predict_cumulative_hazard_function = chf_function,
#'                      predict_function = risk_function)
#'
#' @export
risk_from_chf <- function(predict_cumulative_hazard_function, times) {
    function(model, newdata) rowSums(predict_cumulative_hazard_function(model, newdata, times))
}

#' @keywords internal
add_rug_to_plot <- function(base_plot, rug_df, rug, rug_colors){
    if (rug == "all"){
        return_plot <- with(rug_df, { base_plot +
                geom_rug(data = rug_df[rug_df$statuses == 1,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[1]) +
                geom_rug(data = rug_df[rug_df$statuses == 0,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[2]) })
    } else if (rug == "events") {
        return_plot <- with(rug_df, { base_plot +
                geom_rug(data = rug_df[rug_df$statuses == 1,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[1]) })
    } else if (rug == "censors") {
        return_plot <- with(rug_df, { base_plot +
                geom_rug(data = rug_df[rug_df$statuses == 0,], mapping = aes(x=times, color = statuses), inherit.aes=F, color = rug_colors[2]) })
    } else {
        return_plot <- base_plot
    }
}


#' @keywords internal
calculate_integral <- function(values, times, normalization = "t_max", ...){
    n <- length(values)

    if (is.null(normalization)){
        tmp <- (values[1:(n - 1)] + values[2:n]) * diff(times) / 2
        integrated_metric <- sum(tmp) / (max(times) - min(times))
        return(integrated_metric)
    }
    else if (normalization == "t_max") {
        tmp <- (values[1:(n - 1)] + values[2:n]) * diff(times) / 2
        integrated_metric <- sum(tmp)
        return(integrated_metric/max(times))
    } else if (normalization == "survival"){
        y_true <- list(...)$y_true
        km <- survival::survfit(y_true ~ 1)
        estimator <- stepfun(km$time, c(1, km$surv))

        dwt <- 1 - estimator(times)

        tmp <- (values[1:(n - 1)] + values[2:n]) * diff(dwt) / 2
        integrated_metric <- sum(tmp)
        return(integrated_metric/(1 - estimator(max(times))))
    }
}

# based on iml::order_levels
#' @importFrom stats ecdf xtabs cmdscale
#' @keywords internal
order_levels <- function(data, variable_values, variable_name) {
    feature <- droplevels(variable_values)
    x.count <- as.numeric(table(feature))
    x.prob <- x.count / sum(x.count)
    K <- nlevels(feature)

    dists <- lapply(setdiff(colnames(data), variable_name), function(x) {
        feature.x <- data[, x]
        dists <- expand.grid(levels(feature), levels(feature))
        colnames(dists) <- c("from.level", "to.level")
        if (inherits(feature.x, "factor")) {
            A <- table(feature, feature.x) / x.count
            dists$dist <- rowSums(abs(A[dists[, "from.level"], ] - A[dists[, "to.level"], ])) / 2
        } else {
            quants <- quantile(feature.x, probs = seq(0, 1, length.out = 100), na.rm = TRUE, names = FALSE)
            ecdfs <- data.frame(lapply(levels(feature), function(lev) {
                x.ecdf <- ecdf(feature.x[feature == lev])(quants)
            }))
            colnames(ecdfs) <- levels(feature)
            ecdf.dists.all <- abs(ecdfs[, dists$from.level] - ecdfs[, dists$to.level])
            dists$dist <- apply(ecdf.dists.all, 2, max)
        }
        dists
    })

    dists.cumulated.long <- Reduce(function(d1, d2) {
        d1$dist <- d1$dist + d2$dist
        d1
    }, dists)
    dists.cumulated <- xtabs(dist ~ from.level + to.level, dists.cumulated.long)
    scaled <- cmdscale(dists.cumulated, k = 1)
    order(scaled)
}

