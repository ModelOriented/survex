# survex (development)

* improved and unified API documentation ([#2](https://github.com/ModelOriented/survex/issues/2))
* added references to used methods ([#5](https://github.com/ModelOriented/survex/issues/5))
* changed the package used to draw complex plots from `gridExtra` to `patchwork` ([#7](https://github.com/ModelOriented/survex/pull/7))
* fixed subtitles in plots ([#11](https://github.com/ModelOriented/survex/issues/11))
* fixed calculating of roc curves for classification problems
* added wrapper function for measures provided by `mlr3proba` ([#10](https://github.com/ModelOriented/survex/issues/10))
* created vignette showing how to use `mlr3proba` with `survex`
* fixed incompatibility with new ggplot2 version 3.4 


# survex 0.1.1
* The `survex` package is now public
* `model_parts`, `model_profile`, `predict_parts`, `predict_profile` explanations implemented
* C/D AUC, Brier score and (Harrell's) concordance index performance measures implemented
* Explain methods for `survival`, `ranger`, `randomForestSRC`, `censored` and `mlr3proba` packages.
