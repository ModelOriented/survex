% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_surv_lime.R
\name{plot.surv_lime}
\alias{plot.surv_lime}
\title{Plot SurvLIME Explanations for Survival Models}
\usage{
\method{plot}{surv_lime}(
  x,
  type = "local_importance",
  show_survival_function = TRUE,
  ...,
  title = "SurvLIME",
  subtitle = "default",
  max_vars = 7,
  colors = NULL
)
}
\arguments{
\item{x}{an object of class \code{"surv_lime"} to be plotted}

\item{type}{character, either "coefficients" or "local_importance" (default), selects the type of plot}

\item{show_survival_function}{logical, if the survival function of the explanations should be plotted next to the barplot}

\item{...}{other parameters currently ignored}

\item{title}{character, title of the plot}

\item{subtitle}{character, subtitle of the plot, \code{'default'} automatically generates "created for XXX, YYY models", where XXX and YYY are the explainer labels}

\item{max_vars}{maximum number of variables to be plotted (least important variables are ignored)}

\item{colors}{character vector containing the colors to be used for plotting variables (containing either hex codes "#FF69B4", or names "blue")}
}
\value{
An object of the class \code{ggplot}.
}
\description{
This functions plots objects of class \code{surv_lime} - LIME explanations of survival models
created using \code{predict_parts(..., type="survlime")} function.
}
\examples{
\donttest{
library(survival)
library(survex)

model <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)
exp <- explain(model)

p_parts_lime <- predict_parts(exp, veteran[1, -c(3, 4)], type = "survlime")
plot(p_parts_lime)
}
}
\seealso{
Other functions for plotting 'predict_parts_survival' objects: 
\code{\link{plot.predict_parts_survival}()},
\code{\link{plot.surv_shap}()}
}
\concept{functions for plotting 'predict_parts_survival' objects}
