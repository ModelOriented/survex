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

    if (is.null(colors)) return(colors_discrete_drwhy(n))
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


#' Order levels of a categorical features
#'
#' @description
#' Orders the levels by their similarity in other features. Computes per feature
#' the distance, sums up all distances and does multi-dimensional scaling
#'
#' @details
#' Goal: Compute the distances between two categories.
#' Input: Instances from category 1 and 2
#'
#' 1. For all features, do (excluding the categorical feature for which we are computing the order):
#'  - If the feature is numerical: Take instances from category 1, calculate the
#'  empirical cumulative probability distribution function (ecdf) of the
#'  feature. The ecdf is a function that tells us for a given feature value, how
#'  many values are smaller. Do the same for category 2. The distance is the
#'  absolute maximum point-wise distance of the two ecdf. Practically, this
#'  value is high when the distribution from one category is strongly shifted
#'  far away from the other. This measure is also known as the
#'  Kolmogorov-Smirnov distance
#'  (\url{https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test}).
#'  - If the feature is categorical: Take instances from category 1 and
#'  calculate a table with the relative frequency of each category of the other
#'  feature. Do the same for instances from category 2. The distance is the sum
#'  of the absolute difference of both relative frequency tables.
#' 2. Sum up the distances over all features
#'
#' This algorithm we run for all pairs of categories.
#' Then we have a k times k matrix, when k is the number of categories, where
#' each entry is the distance between two categories. Still not enough to have a
#' single order, because, a (dis)similarity tells you the pair-wise distances,
#' but does not give you a one-dimensional ordering of the classes. To kind of
#' force this thing into a single dimension, we have to use a dimension
#' reduction trick called multi-dimensional scaling. This can be solved using
#' multi-dimensional scaling, which takes in a distance matrix and returns a
#' distance matrix with reduced dimension. In our case, we only want 1 dimension
#' left, so that we have a single ordering of the categories and can compute the
#' accumulated local effects. After reducing it to a single ordering, we are
#' done and can use this ordering to compute ALE. This is not the Holy Grail how
#' to order the factors, but one possibility.
#'
#' @param data_dt data.frame with the training data
#' @param feature_name the name of the categorical feature
#' @return the order of the levels (not levels itself)
#' @keywords internal
order_levels <- function(data, variable) {
    data[, variable] <- droplevels(data[, variable])
    feature <- data[, variable]
    x.count <- as.numeric(table(data[, variable]))
    x.prob <- x.count / sum(x.count)
    K <- nlevels(data[, variable])

    dists <- lapply(setdiff(colnames(data), variable), function(x) {
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

