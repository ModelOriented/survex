# survex: Explainable Machine Learning in Survival Analysis <img src="man/figures/survex.png" align="right" width="150px"/>

[![R-CMD-check](https://github.com/ModelOriented/survex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ModelOriented/survex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ModelOriented/survex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ModelOriented/survex?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/survex)](https://cran.r-project.org/package=survex)

## Overview 

Survival analysis is a task dealing with time-to-event prediction. Aside from the well understood models, many others have recently emerged, however most of them lack interpretability as they are black-box models. Due to the unusual, functional type of prediction (either in the form of survival function or cumulative hazard function) standard model agnostic explanations cannot be applied directly.

The `survex` package provides model agnostic explanations for survival models. If you're unfamiliar with model agnostic explanations, consider looking at the [Explanatory Model Analysis](https://ema.drwhy.ai/) e-book, as most of the methods included in this package are extensions of those described in the book for models with functional outputs. 

The main function `explain()` creates a standardized wrapper for a model, which is further used for calculating predictions. If you work with models from `mlr3proba`, `censored`, `ranger`, `randomForestSRC` or `survival` packages, creating explainers is automated, most often you only need to supply the `model` parameter to the `explain()` function.

However, an explainer can be created for **any** survival model, using the `explain_survival()` function by passing `model`, `data`, `y`, and `predict_survival_function` arguments.


## Installation

The package is available on [CRAN](https://cran.r-project.org/package=survex):

```r
install.packages("survex")
```

The latest development version can be installed from GitHub using `devtools::install_github()`:

```r
devtools::install_github("https://github.com/ModelOriented/survex")
```


## Usage

[![`survex` usage cheatsheet](man/figures/cheatsheet.png)](https://github.com/ModelOriented/survex/blob/main/misc/cheatsheet.pdf)


## Related work

- H. Ishwaran, U. B. Kogalur, E. H. Blackstone, M. S. Lauer. [Random survival forests](https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F08-AOAS169). Annals of Applied Statistics, 2008.
- M. S. Kovalev, L. V. Utkin, E. M. Kasimov. [SurvLIME: A method for explaining machine learning survival models](https://doi.org/10.1016/j.knosys.2020.106164). Knowledge-Based Systems, 2020.
- R. Sonabend, F. J. Király, A. Bender, B. Bischl, M. Lang. [mlr3proba: an R package for machine learning in survival analysis](https://doi.org/10.1093/bioinformatics/btab039). Bioinformatics, 2021.
- E. Hvitfeldt, H. Frick. [censored: 'parsnip' Engines for Survival Models](https://github.com/tidymodels/censored). CRAN v0.1.0, 2022.
- M. Krzyziński, M. Spytek, H. Baniecki, P. Biecek. [SurvSHAP(t): Time-dependent explanations of machine learning survival models](https://arxiv.org/abs/2208.11080). arXiv preprint arXiv:2208.11080, 2022.
