.survex.env <- new.env()

#' Default Theme for survex plots
#'
#' @param default_theme object - string ("drwhy" or "ema") or an object of ggplot theme class. Will be applied by default by `survex` to all horizontal plots
#' @param default_theme_vertical object - string ("drwhy" or "ema") or an object of ggplot theme class. Will be applied by default by `survex` to all vertical plots
#' @return list with current default themes
#' @examples
#' old <- set_theme_survex("ema")
#' \donttest{
#' library(survival)
#' library(survex)
#'
#' model <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
#' exp <- explain(model)
#'
#' p_parts_lime <- predict_parts(exp, veteran[1, -c(3, 4)], type = "survlime")
#' old <- set_theme_survex("drwhy")
#' plot(p_parts_lime)
#' old <- set_theme_survex(ggplot2::theme_void(), ggplot2::theme_void())
#' plot(p_parts_lime)
#'}
#'
#' @importFrom DALEX theme_drwhy theme_drwhy_vertical theme_ema theme_ema_vertical
#'
#' @export
#' @rdname theme_survex
#'
set_theme_survex <- function(default_theme = "drwhy", default_theme_vertical = default_theme) {
    # it should be either name or theme object
    if (!(any(
        class(default_theme) %in% c("character", "theme")
    )))
        stop("The 'default_theme' shall be either character 'drwhy'/'ema' or ggplot2::theme object")
    if (!(any(
        class(default_theme_vertical) %in% c("character", "theme")
    )))
        stop("The 'default_theme_vertical' shall be either character 'drwhy'/'ema' or ggplot2::theme object")

    # get default themes
    old <- .survex.env$default_themes

    # set themes
    if (is.character(default_theme)) {
        # from name
        switch (default_theme,
                drwhy = {.survex.env$default_themes <- list(default = theme_drwhy(), vertical = theme_drwhy_vertical())},
                ema = {.survex.env$default_themes <- list(default = theme_ema(), vertical = theme_ema_vertical())},
                stop("Only 'drwhy' or 'ema' names are allowed")
        )
    } else {
        # from themes (ggplot2 objects)
        .survex.env$default_themes <- list(default = default_theme, vertical = default_theme_vertical)
    }

    # return old themes
    old
}


#' @export
#' @rdname theme_survex
theme_default_survex <- function() {
    if (!exists("default_themes", envir = .survex.env))
        return(theme_drwhy())

    .survex.env$default_themes[[1]]
}

#' @export
#' @rdname theme_survex
theme_vertical_default_survex <- function() {
    if (!exists("default_themes", envir = .survex.env))
        return(theme_drwhy_vertical())

    .survex.env$default_themes[[2]]
}
