# survex: Explainable Machine Learning in Survival Analysis <img src="man/figures/survex.png" align="right" width="150px"/>

[![R-CMD-check](https://github.com/ModelOriented/survex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ModelOriented/survex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ModelOriented/survex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ModelOriented/survex?branch=main)

## Overview 

Survival analysis is a task dealing with time-to-event prediction. Aside from the well understood models, many others have recently emerged, however most of them lack interpretability as they are black-box models. Due to the unusual, functional type of prediction (either in the form of survival function or cumulative hazard function) standard model agnostic explanations cannot be applied directly.

The `survex` package provides model agnostic explanations for survival models. If you're unfamiliar with model agnostic explanations, consider looking at the [Explanatory Model Analysis](https://ema.drwhy.ai/) e-book, as most of the methods included in this package are extensions of those described in the book for models with functional outputs. 

The main function `explain()` creates a standardized wrapper for a model, which is further used for calculating predictions. If you work with models from `mlr3proba`, `censored`, `ranger`, `randomForestSRC` or `survival` packages, creating explainers is automated, most often you only need to supply the `model` parameter to the `explain()` function.

However, an explainer can be created for **any** survival model, using the `explain_survival()` function by passing `model`, `data`, `y`, and `predict_survival_function` arguments.

## Installation

The package can be installed from github using `devtools::install_github()`:

```
devtools::install_github("https://github.com/ModelOriented/survex")
```

## Usage

![`survex` usage cheatsheet](man/figures/cheatsheet.png)

## Related work

- mlr3proba: an R package for machine learning in survival analysis [[Bioinformatics]](https://academic.oup.com/bioinformatics/article/37/17/2789/6125361)
- censored: an extention package for parsnip containing survival analysis models [[censored]](https://censored.tidymodels.org/) 
