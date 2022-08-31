#' Extract info from model
#'
#' This generic function let user extract base information about model. The function returns a named list of class \code{model_info} that
#' contain information about package of model, version and task type. For wrappers like \code{mlr} or \code{parsnip} both, package and wrapper
#' information are stored
#'
#' @param model - model object
#' @param ... - another arguments
#'
#' @details
#' Currently supported packages are:
#' \itemize{
#' \item class \code{coxph} - Cox proportional hazards regression model created with \pkg{survival} package
#' \item class \code{model_fit} - models created with \pkg{parsnip} package
#' \item class \code{ranger} - random survival forest models created with \pkg{ranger} package
#' \item class \code{rfsrc} - random forest models created with \pkg{randomForestSRC} package
#' }
#'
#' @return A named list of class \code{model_info}
#'
#' @rdname surv_model_info
#' @export
#'
#' @examples
#' library(survival)
#' library(survex)
#' cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran,
#'                        model = TRUE, x = TRUE, y = TRUE)
#' surv_model_info(cph)
#'
#' \donttest{
#' library(ranger)
#' rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran,
#'  num.trees = 50, mtry = 3, max.depth = 5)
#' surv_model_info(rsf_ranger)
#' }
#'
surv_model_info <- function(model, ...)
    UseMethod("surv_model_info")


#' @rdname surv_model_info
#' @export
surv_model_info.coxph <- function(model, ...) {
    type <- "survival"
    package <- "survival"
    ver <- get_pkg_ver_safe(package)
    model_info <- list(package = package, ver = ver, type = type)
    class(model_info) <- "model_info"
    model_info
}

#' @rdname surv_model_info
#' @export
surv_model_info.rfsrc <- function(model, ...) {
    type <- "survival"
    package <- "randomForestSRC"
    ver <- get_pkg_ver_safe(package)
    model_info <- list(package = package, ver = ver, type = type)
    class(model_info) <- "model_info"
    model_info
}

#' @rdname surv_model_info
#' @export
surv_model_info.ranger <- function(model, ...) {
    type <- "survival"
    package <- "ranger"
    ver <- get_pkg_ver_safe(package)
    model_info <- list(package = package, ver = ver, type = type)
    class(model_info) <- "model_info"
    model_info
}

#' @rdname surv_model_info
#' @export
surv_model_info.model_fit <- function(model, ...) {
    type <- "survival"
    package_wrapper <- "parsnip"
    ver_wrapper <- get_pkg_ver_safe(package_wrapper)
    package <- model$spec$method$libs
    ver <- get_pkg_ver_safe(package)
    model_info <- list(package = c(wrapper = package_wrapper, package = package), ver = c(wrapper = ver_wrapper, package = ver), type = type)
    class(model_info) <- "model_info"
    model_info
}

#' @rdname surv_model_info
#' @export
surv_model_info.cph <- function(model, ...) {
    type <- "survival"
    package_wrapper <- "rms"
    ver_wrapper <- get_pkg_ver_safe(package_wrapper)
    package <- "survival"
    ver <- get_pkg_ver_safe(package)
    model_info <- list(package = c(wrapper = package_wrapper, package = package), ver = c(wrapper = ver_wrapper, package = ver), type = type)
    class(model_info) <- "model_info"
    model_info
}


#' @rdname surv_model_info
#' @export
surv_model_info.LearnerSurv <- function(model, ...) {
    type <- "survival"
    package <- "mlr3proba"
    ver <- get_pkg_ver_safe(package)
    model_info <- list(package = package, ver = ver, type = type)
    class(model_info) <- "model_info"
    model_info
}



#' @rdname surv_model_info
#' @export
surv_model_info.default <- function(model, ...) {
    type <- "survival"
    package <- paste("Model of class:", class(model), "package unrecognized")
    ver <- "Unknown"
    model_info <- list(package = package, ver = ver, type = type)
    class(model_info) <- "model_info"
    model_info
}

get_pkg_ver_safe <- function(package) {
    ver <- try(as.character(utils::packageVersion(package)), silent = TRUE)
    if (inherits(ver, "try-error")) {
        ver <- "Unknown"
    }
    ver
}
