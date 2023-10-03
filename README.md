# survex: Explainable Machine Learning in Survival Analysis <img src="man/figures/survex.png" align="right" width="150px"/>

[![R-CMD-check](https://github.com/ModelOriented/survex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ModelOriented/survex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ModelOriented/survex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ModelOriented/survex?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/survex)](https://cran.r-project.org/package=survex)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/survex?color=orange)](https://cranlogs.r-pkg.org/badges/grand-total/survex)
[![DrWhy-BackBone](https://img.shields.io/badge/DrWhy-BackBone-373589)](http://drwhy.ai/#BackBone)


## Overview 

Survival analysis is a task dealing with time-to-event prediction. Aside from the well-understood models like CPH, many more complex models have recently emerged, but most lack interpretability. Due to a functional type of prediction, either in the form of survival function or cumulative hazard function, standard model-agnostic explanations cannot be applied directly.

The `survex` package provides model-agnostic explanations for machine learning survival models. It is based on the [`DALEX` package](https://github.com/ModelOriented/DALEX). If you're unfamiliar with explainable machine learning, consider referring to the [Explanatory Model Analysis](https://ema.drwhy.ai) book -- most of the methods included in `survex` extend these described in EMA and implemented in `DALEX` but to models with functional output. 

The main `explain()` function uses a model and data to create a standardized `explainer` object, which is further used as an interface for calculating predictions. We automate creating explainers from the following packages: `mlr3proba`, `censored`, `ranger`, `randomForestSRC`, and `survival`. **Raise an Issue on GitHub if you find models from other packages that we can incorporate into the `explain()` interface.**

Note that an explainer can be created for **any** survival model, using the `explain_survival()` function by passing `model`, `data`, `y`, and `predict_survival_function` arguments.


## Installation

The package is available on [CRAN](https://cran.r-project.org/package=survex):

```r
install.packages("survex")
```

The latest development version can be installed from GitHub using `devtools::install_github()`:

```r
devtools::install_github("https://github.com/ModelOriented/survex")
```

## Simple demo

```r
library("survex")
library("survival")
library("ranger")

# create a model
model <- ranger(Surv(time, status) ~ ., data = veteran)

# create an explainer
explainer <- explain(model, 
                     data = veteran[, -c(3, 4)],
                     y = Surv(veteran$time, veteran$status))

# evaluate the model
model_performance(explainer)

# visualize permutation-based feature importance
plot(model_parts(explainer))

# explain one prediction with SurvSHAP(t)
plot(predict_parts(explainer, veteran[1, -c(3, 4)]))
```

## Functionalities and roadmap

Existing functionalities:
- [x] unified prediction interface using the explainer object - `predict()`
- [x] calculation of performance metrics (Brier Score, Time-dependent C/D AUC, metrics from `mlr3proba`) - `model_performance()`
- [x] calculation of feature importance (Permutation Feature Importance - PFI) - `model_parts()`
- [x] calculation of partial dependence (Partial Dependence Profiles - PDP, Accumulated Local Effects - ALE) - `model_profile()`
- [x] calculation of 2-dimensional partial dependence (2D PDP, 2D ALE) - `model_profile_2d()`
- [x] calculation of local feature attributions (SurvSHAP(t), SurvLIME) - `predict_parts()`
- [x] calculation of local ceteris paribus explanations (Ceteris Paribus profiles - CP/ Individual Conditional Expectations - ICE) - `predict_profile()`
- [x] calculation of global feature attributions using SurvSHAP(t) - `model_survshap()`

Currently in develompment:
- [ ] ...

Future plans:
- [ ] ... (raise an Issue on GitHub if you have any suggestions)

## Usage

[![`survex` usage cheatsheet](man/figures/cheatsheet.png)](https://github.com/ModelOriented/survex/blob/main/misc/cheatsheet.pdf)


## Citation

If you use `survex`, please cite [our preprint](https://arxiv.org/abs/2308.16113):

> M. Spytek, M. Krzyziński, S. H. Langbein, H. Baniecki, M. N. Wright, P. Biecek. *survex: an R package for explaining machine learning survival models*. **arXiv preprint arXiv:2308.16113**, 2023.

```
@article{spytek2023survex,
    title   = {{survex: an R package for explaining machine learning survival models}},
    author  = {Mikołaj Spytek and Mateusz Krzyziński and Sophie Hanna Langbein and
               Hubert Baniecki and Marvin N. Wright and Przemysław Biecek},
    journal = {arXiv preprint arXiv:2308.16113},
    year    = {2023}
}
```

## Applications of `survex`

- H. Baniecki, B. Sobieski, P. Bombiński, P. Szatkowski, P. Biecek. [Hospital Length of Stay Prediction Based on Multi-modal Data towards Trustworthy Human-AI Collaboration in Radiomics](https://arxiv.org/abs/2303.09817). *International Conference on Artificial Intelligence in Medicine*, 2023.
- W. Chen, B. Zhou, C. Y. Jeon, F. Xie, Y-C. Lin, R. K. Butler, Y. Zhou, T. Q. Luong, E. Lustigova, J. R. Pisegna, B. U. Wu. [Machine learning versus regression for prediction of sporadic pancreatic cancer](https://doi.org/10.1016/j.pan.2023.04.009). *Pancreatology*, 2023.
- M. Nachit, Y. Horsmans, R. M. Summers, I. A. Leclercq, P. J. Pickhardt. [AI-based CT Body Composition Identifies Myosteatosis as Key Mortality Predictor in Asymptomatic Adults](https://doi.org/10.1148/radiol.222008). *Radiology*, 2023.
- R. Passera, S. Zompi, J. Gill, A. Busca. [Explainable Machine Learning (XAI) for Survival in Bone Marrow Transplantation Trials: A Technical Report](https://doi.org/10.3390/biomedinformatics3030048). *BioMedInformatics*, 2023.
- P. Donizy, M. Spytek, M. Krzyziński, K. Kotowski, A. Markiewicz, B. Romanowska-Dixon, P. Biecek, M. P. Hoang. [Ki67 is a better marker than PRAME in risk stratification of BAP1-positive and BAP1-loss uveal melanomas](http://dx.doi.org/10.1136/bjo-2023-323816). *British Journal of Ophthalmology*, 2023.
- X. Qi, Y. Ge, A. Yang, Y. Liu, Q. Wang & G. Wu. [Potential value of mitochondrial regulatory pathways in the clinical application of clear cell renal cell carcinoma: a machine learning-based study](https://doi.org/10.1007/s00432-023-05393-8). *Journal of Cancer Research and Clinical Oncology*, 2023.

- Share it with us!

## Related work

- H. Ishwaran, U. B. Kogalur, E. H. Blackstone, M. S. Lauer. [Random survival forests](https://doi.org/10.1214/08-AOAS169). *Annals of Applied Statistics*, 2008.
- A. Grudziąż, A. Gosiewska, P. Biecek. [survxai: an R package for structure-agnostic explanations of survival models](https://doi.org/10.21105/joss.00961). *Journal of Open Source Software*, 2018.
- M. S. Kovalev, L. V. Utkin, E. M. Kasimov. [SurvLIME: A method for explaining machine learning survival models](https://doi.org/10.1016/j.knosys.2020.106164). *Knowledge-Based Systems*, 2020.
- R. Sonabend, F. J. Király, A. Bender, B. Bischl, M. Lang. [mlr3proba: an R package for machine learning in survival analysis](https://doi.org/10.1093/bioinformatics/btab039). *Bioinformatics*, 2021.
- E. Hvitfeldt, H. Frick. [censored: 'parsnip' Engines for Survival Models](https://github.com/tidymodels/censored). *CRAN v0.1.0*, 2022.
- M. Krzyziński, M. Spytek, H. Baniecki, P. Biecek. [SurvSHAP(t): Time-dependent explanations of machine learning survival models](https://doi.org/10.1016/j.knosys.2022.110234). *Knowledge-Based Systems*, 2023.
