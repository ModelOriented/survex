% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_model_parts_survival.R
\name{plot.model_parts_survival}
\alias{plot.model_parts_survival}
\title{Plot Model Parts for Survival Models}
\usage{
\method{plot}{model_parts_survival}(x, ...)
}
\arguments{
\item{x}{an object of class \code{"model_parts_survival"} to be plotted}

\item{...}{additional parameters passed to the \code{plot.surv_feature_importance} function}
}
\value{
An object of the class \code{ggplot}.
}
\description{
This function is a wrapper for plotting \code{model_parts} objects created for survival models
using the \code{model_parts()} function.
}
\section{Plot options}{

\itemize{
\item \code{title} - character, title of the plot
\item \code{subtitle} - character, subtitle of the plot, if \code{NULL} automatically generated as "created for XXX, YYY models", where XXX and YYY are explainer labels
\item \code{max_vars} - maximum number of variables to be plotted (least important variables are ignored)
\item \code{colors} - character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")
\item \code{rug} - character, one of \code{"all"}, \code{"events"}, \code{"censors"}, \code{"none"} or \code{NULL}. Which times to mark on the x axis in \code{geom_rug()}.
\item \code{rug_colors} - character vector containing two colors (containing either hex codes "#FF69B4", or names "blue"). The first color (red by default) will be used to mark event times, whereas the second (grey by default) will be used to mark censor times.
}
}

\examples{
\donttest{
library(survival)
library(survex)

model <- coxph(Surv(time, status) ~ ., data = veteran, x = TRUE, model = TRUE, y = TRUE)
explainer <- explain(model)

mp <- model_parts(explainer)

plot(mp)
}
}
\seealso{
Other functions for plotting 'model_parts_survival' objects: 
\code{\link{plot.surv_feature_importance}()}
}
\concept{functions for plotting 'model_parts_survival' objects}
