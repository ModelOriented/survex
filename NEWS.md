# survex (development)

* change behaviour of `categorical_variables` argument in `model_parts()` and `predict_parts()`. If it containes variable names not present in the `variables` argument, they will be added at the end. ([#39](https://github.com/ModelOriented/survex/issues/39))
* fix invalid color palette order in plot feature importance
* fix predict_parts survshap running out of memory with more than 16 variables ([#25](https://github.com/ModelOriented/survex/issues/25))

# survex 0.2.2

* improved and unified API documentation ([#2](https://github.com/ModelOriented/survex/issues/2))
* added references to used methods ([#5](https://github.com/ModelOriented/survex/issues/5))
* changed the package used to draw complex plots from `gridExtra` to `patchwork` ([#7](https://github.com/ModelOriented/survex/pull/7))
* fixed subtitles in plots ([#11](https://github.com/ModelOriented/survex/issues/11))
* fixed calculating of ROC curves for classification problems
([#17](https://github.com/ModelOriented/survex/issues/17))
* added wrapper function for measures provided by `mlr3proba` ([#10](https://github.com/ModelOriented/survex/issues/10))
* created vignette showing how to use `mlr3proba` with `survex`
* fixed incompatibility with new ggplot2 version 3.4 
* added function for creating integrated versions of time-dependent metrics ([#9](https://github.com/ModelOriented/survex/issues/9))
* move `ingredients` from imports to suggests


# survex 0.1.1
* The `survex` package is now public
* `model_parts`, `model_profile`, `predict_parts`, `predict_profile` explanations implemented
* C/D AUC, Brier score and (Harrell's) concordance index performance measures implemented
* Explain methods for `survival`, `ranger`, `randomForestSRC`, `censored` and `mlr3proba` packages.
