# survex (development)

* *breaking change:* refactored the structure of `model_performance_survival` object - calculated metrics are now in a `$result` list.
* fixed invalid color palette order in plot feature importance
* fixed predict_parts survshap running out of memory with more than 16 variables ([#25](https://github.com/ModelOriented/survex/issues/25))
* added `max_vars` parameter for predict_parts explanations ([#27](https://github.com/ModelOriented/survex/issues/27))
* set `max_vars` to 7 for every method  
* refactored survshap code ([#29](https://github.com/ModelOriented/survex/issues/29), [#30](https://github.com/ModelOriented/survex/issues/30), [#43](https://github.com/ModelOriented/survex/issues/43))
* fixed survshap error when target columns named different than time and status ([#44](https://github.com/ModelOriented/survex/issues/44))
* fixed survlime error when all variables are categorical ([#46](https://github.com/ModelOriented/survex/issues/46))
* fixed subtitles in feature importance plots ([#11](https://github.com/ModelOriented/survex/issues/11))
* added the possibility to set themes with `set_theme_survex()` ([#32](https://github.com/ModelOriented/survex/issues/32))
* added the possibility of plotting multiple `predict_parts()` and `model_parts()` explanations in one plot ([#12](https://github.com/ModelOriented/survex/issues/12))
* fixed the x axis of plots (it now starts from 0) ([#37](https://github.com/ModelOriented/survex/issues/37))
* added geom_rug() to all time-dependent plots, marking event and censoring times ([#35](https://github.com/ModelOriented/survex/issues/35))
* refactored `surv_feature_importance.R` - change auxiliary columns to include `_` in their name. Necessary changes also done to plotting and printing functions. ([#28](https://github.com/ModelOriented/survex/issues/28))
* changed default `type` argument of `model_parts()` to `"difference"` ([#33](https://github.com/ModelOriented/survex/issues/33))
* refactored integration of metrics ([#31](https://github.com/ModelOriented/survex/issues/31))
* changed behaviour of `categorical_variables` argument in `model_parts()` and `predict_parts()`. If it contains variable names not present in the `variables` argument, they will be added at the end. ([#39](https://github.com/ModelOriented/survex/issues/39))
* added ROC AUC calculation and plotting for selected timepoints in `model_performance()` ([#22](https://github.com/ModelOriented/survex/issues/22))


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
