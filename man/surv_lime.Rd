% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surv_lime.R
\name{surv_lime}
\alias{surv_lime}
\title{Helper functions for \code{predict_parts.R}}
\usage{
surv_lime(
  explainer,
  new_observation,
  ...,
  N = 100,
  distance_metric = "euclidean",
  kernel_width = NULL,
  sampling_method = "gaussian",
  sample_around_instance = TRUE,
  max_iter = 10000,
  categorical_variables = NULL,
  k = 1 + 1e-04
)
}
\arguments{
\item{explainer}{an explainer object - model preprocessed by the \code{explain()} function}

\item{new_observation}{a new observation for which predictions need to be explained}

\item{...}{additional parameters, passed to internal functions}

\item{N}{a positive integer, number of observations generated in the neighbourhood}

\item{distance_metric}{character, name of the distance metric to be used, only \code{"euclidean"} is implemented}

\item{kernel_width}{a numeric or \code{"silverman"}, parameter used for calculating weights, by default it's \code{sqrt(ncol(data)*0.75)}. If \code{"silverman"} the kernel width is calculated using the method proposed by Silverman and used in the SurvLIMEpy Python package.}

\item{sampling_method}{character, name of the method of generating neighbourhood, only \code{"gaussian"} is implemented}

\item{sample_around_instance}{logical, if the neighbourhood should be generated with the new observation as the center (default), or should the mean of the whole dataset be used as the center}

\item{max_iter}{a numeric, maximal number of iteration for the optimization problem}

\item{categorical_variables}{character vector, names of variables that should be treated as categories (factors are included by default)}

\item{k}{a small positive number > 1, added to chf before taking log, so that weigths aren't negative}
}
\value{
A list, with the SurvLIME result in the \verb{$result} field.
}
\description{
Helper functions for \code{predict_parts.R}
}
\section{References}{

\itemize{
\item [1] Kovalev, Maxim S., et al. \href{https://www.sciencedirect.com/science/article/pii/S0950705120304044?casa_token=6e9cyk_ji3AAAAAA:tbqo33MsZvNC9nrSGabZdLfPtZTsvsvZTHYQCM2aEhumLI5D46U7ovhr37EaYUhmKZrw45JzDhg}{"SurvLIME: A method for explaining machine learning survival models."} Knowledge-Based Systems 203 (2020): 106164.
\item [2] Pachón-García, Cristian, et al. \href{https://www.sciencedirect.com/science/article/pii/S095741742302122X}{"SurvLIMEpy: A Python package implementing SurvLIME."} Expert Systems with Applications 237 (2024): 121620.
}
}

\keyword{internal}
