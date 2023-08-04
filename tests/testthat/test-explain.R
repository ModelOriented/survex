test_that("coxph prediction functions work correctly", {
    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL
    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )


    coxph_explainer <-
        explain(cox_rotterdam_rec,
                y = survival::Surv(rotterdam$rtime, rotterdam$recur),
                verbose = FALSE)

    times <- coxph_explainer$times


    expect_equal(
        coxph_explainer$predict_survival_function(cox_rotterdam_rec, rotterdam[c(1, 2, 3), ], times),
        pec::predictSurvProb(cox_rotterdam_rec, rotterdam[c(1, 2, 3), ], times)
    )
    expect_equal(
        coxph_explainer$predict_survival_function(cox_rotterdam_rec, rotterdam[4, ], times),
        pec::predictSurvProb(cox_rotterdam_rec, rotterdam[4, ], times)
    )

    expect_equal(
        coxph_explainer$predict_cumulative_hazard_function(cox_rotterdam_rec, rotterdam[c(1, 2, 3), ], times),
        survival_to_cumulative_hazard(pec::predictSurvProb(cox_rotterdam_rec, rotterdam[c(1, 2, 3), ], times))
    )
    expect_equal(
        coxph_explainer$predict_cumulative_hazard_function(cox_rotterdam_rec, rotterdam[4, ], times),
        survival_to_cumulative_hazard(pec::predictSurvProb(cox_rotterdam_rec, rotterdam[4, ], times))
    )

    expect_equal(
        coxph_explainer$predict_function(cox_rotterdam_rec, rotterdam[c(1, 2, 3), ]),
        predict(cox_rotterdam_rec, rotterdam[c(1, 2, 3), ], type = "risk")
    )
    expect_equal(
        coxph_explainer$predict_function(cox_rotterdam_rec, rotterdam[4, ]),
        predict(cox_rotterdam_rec, rotterdam[4, ], type = "risk")
    )
})

test_that("ranger prediction functions work correctly", {
    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL
    rsf_rotterdam_rec <-
        ranger::ranger(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            respect.unordered.factors = TRUE,
            num.trees = 100,
            mtry = 3,
            max.depth = 4
        )

    ranger_explainer <-
        explain(rsf_rotterdam_rec,
                data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death", "rtime", "recur")],
                y = survival::Surv(rotterdam$rtime, rotterdam$recur),
                verbose = FALSE
                )

    times <- ranger_explainer$times


    expect_equal(
        ranger_explainer$predict_survival_function(rsf_rotterdam_rec, rotterdam[c(1, 2, 3), ], times),
        {
            prediction <- predict(rsf_rotterdam_rec, rotterdam[c(1, 2, 3), ])
            return_matrix <- matrix(nrow = 3, ncol = length(times))

            for (i in 1:3) {
                stepfunction <- stepfun(prediction$unique.death.times, c(1, prediction$survival[i, ]))
                return_matrix[i, ] <- stepfunction(times)
            }
            return_matrix
        }
    )
    expect_equal(
        ranger_explainer$predict_survival_function(rsf_rotterdam_rec, rotterdam[4, ], times),
        {
            prediction <- predict(rsf_rotterdam_rec, rotterdam[4, ])
            return_matrix <- matrix(nrow = 1, ncol = length(times))
            stepfunction <- stepfun(prediction$unique.death.times, c(1, prediction$survival))
            return_matrix[1, ] <- stepfunction(times)
            return_matrix

        }
    )

    expect_equal(
        ranger_explainer$predict_cumulative_hazard_function(rsf_rotterdam_rec, rotterdam[c(1, 2, 3), ], times),
        {
            prediction <- predict(rsf_rotterdam_rec, rotterdam[c(1, 2, 3), ])
            return_matrix <- matrix(nrow = 3, ncol = length(times))

            for (i in 1:3) {
                stepfunction <- stepfun(prediction$unique.death.times, c(1, prediction$chf[i, ]))
                return_matrix[i, ] <- stepfunction(times)
            }

            return_matrix
        }
    )
    expect_equal(
        ranger_explainer$predict_cumulative_hazard_function(rsf_rotterdam_rec, rotterdam[4, ], times),
        {
            prediction <- predict(rsf_rotterdam_rec, rotterdam[4, ])
            return_matrix <- matrix(nrow = 1, ncol = length(times))

            stepfunction <- stepfun(prediction$unique.death.times, c(1, prediction$chf))
            return_matrix[1, ] <- stepfunction(times)


            return_matrix
        }
    )

    expect_equal(ranger_explainer$predict_function(rsf_rotterdam_rec, rotterdam[c(1, 2, 3), ]),
                 {
                     prediction <- predict(rsf_rotterdam_rec, rotterdam[c(1, 2, 3), ])
                     return_matrix <- matrix(nrow = 3, ncol = length(times))

                     for (i in 1:3) {
                         stepfunction <-
                             stepfun(prediction$unique.death.times, c(1, prediction$chf[i, ]))
                         return_matrix[i, ] <- stepfunction(times)
                     }

                     rowSums(return_matrix)
                 })
    expect_equal(ranger_explainer$predict_function(rsf_rotterdam_rec, rotterdam[4, ]),
                 {
                     prediction <- predict(rsf_rotterdam_rec, rotterdam[4, ])
                     return_matrix <- matrix(nrow = 1, ncol = length(times))

                     stepfunction <- stepfun(prediction$unique.death.times, c(1, prediction$chf))
                     return_matrix[1, ]  <- stepfunction(times)
                     rowSums(return_matrix)
                 })


})

test_that("rsfrc prediction functions work correctly", {
    colon <- survival::colon
    colon <- colon[colon$etype == 1, ]
    colon$study <- NULL
    colon$id <- NULL
    colon$extent <- as.factor(colon$extent)
    colon <- na.omit(colon)

    rsf_colon_rec <-
        randomForestSRC::rfsrc(
            Surv(time, status) ~ rx + sex + age + obstruct + perfor +
                adhere + nodes + differ +
                extent + surg + node4,
            data = colon,
            subset = colon$etype == 1
        )

    rsf_explainer <-
        explain(rsf_colon_rec, y = survival::Surv(colon$time, colon$status), verbose = FALSE)

    times <- rsf_explainer$times


    expect_equal(rsf_explainer$predict_survival_function(rsf_colon_rec, colon[c(11, 12, 13), ], times),
                 {
                     prediction <- predict(rsf_colon_rec, colon[c(11, 12, 13), ])
                     return_matrix <- matrix(nrow = 3, ncol = length(times))

                     for (i in 1:3) {
                         stepfunction <- stepfun(prediction$time.interest, c(1, prediction$survival[i, ]))
                         return_matrix[i, ] <- stepfunction(times)
                     }
                     return_matrix
                 })

    expect_equal(rsf_explainer$predict_survival_function(rsf_colon_rec, colon[14, ], times),
                 {
                     prediction <- predict(rsf_colon_rec, colon[14, ])
                     return_matrix <- matrix(nrow = 1, ncol = length(times))

                     stepfunction <- stepfun(prediction$time.interest, c(1, prediction$survival))
                     return_matrix[1, ] <- stepfunction(times)

                     return_matrix
                 })

    expect_equal(rsf_explainer$predict_cumulative_hazard_function(rsf_colon_rec, colon[c(11, 12, 13), ], times),
                 {
                     prediction <- predict(rsf_colon_rec, colon[c(11, 12, 13), ])
                     return_matrix <- matrix(nrow = 3, ncol = length(times))

                     for (i in 1:3) {
                         stepfunction <- stepfun(prediction$time.interest, c(1, prediction$chf[i, ]))
                         return_matrix[i, ] <- stepfunction(times)
                     }
                     return_matrix
                 })

    expect_equal(rsf_explainer$predict_cumulative_hazard_function(rsf_colon_rec, colon[14, ], times),
                 {
                     prediction <- predict(rsf_colon_rec, colon[14, ])
                     return_matrix <- matrix(nrow = 1, ncol = length(times))

                     stepfunction <- stepfun(prediction$time.interest, c(1, prediction$chf))
                     return_matrix[1, ] <- stepfunction(times)

                     return_matrix
                 })

    expect_equal(rsf_explainer$predict_function(rsf_colon_rec, colon[c(11, 12, 13), ]),
                 {
                     prediction <- predict(rsf_colon_rec, colon[c(11, 12, 13), ])
                     return_matrix <- matrix(nrow = 3, ncol = length(times))

                     for (i in 1:3) {
                         stepfunction <- stepfun(prediction$time.interest, c(1, prediction$chf[i, ]))
                         return_matrix[i, ] <- stepfunction(times)
                     }
                     rowSums(return_matrix)
                 })
    expect_equal(rsf_explainer$predict_function(rsf_colon_rec, colon[14, ]),
                 {
                     prediction <- predict(rsf_colon_rec, colon[14, ])
                     return_matrix <- matrix(nrow = 1, ncol = length(times))

                     stepfunction <- stepfun(prediction$time.interest, c(1, prediction$chf))
                     return_matrix[1, ] <- stepfunction(times)

                     rowSums(return_matrix)
                 })
})


test_that("automated `y` and `data` sourcing works", {

    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL

    colon <- survival::colon
    colon$study <- NULL
    colon$id <- NULL
    colon$extent <- as.factor(colon$extent)
    colon <- na.omit(colon)

    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )

    rsf_rotterdam_rec <- ranger::ranger(survival::Surv(rtime, recur) ~ .,
                                        data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
                                        respect.unordered.factors = TRUE,
                                        num.trees = 100,
                                        mtry = 3,
                                        max.depth = 4)

    rsf_colon_rec <- randomForestSRC::rfsrc(Surv(time, status) ~ rx + sex + age + obstruct + perfor +
                                                adhere + nodes + differ + extent + surg + node4,
                                            data = colon, subset = colon$etype == 1)

    coxph_explainer <- explain(cox_rotterdam_rec, verbose = FALSE)
    y <- coxph_explainer$y
    dimnames(y)[[1]] <- list()

    expect_equal(y, survival::Surv(rotterdam$rtime, rotterdam$recur), ignore_attr = TRUE)
    expect_equal(predict(cox_rotterdam_rec, coxph_explainer$data), predict(cox_rotterdam_rec, rotterdam))

    rsf_explainer <- explain(rsf_colon_rec, verbose = FALSE)

    expect_equal(rsf_explainer$y, survival::Surv(rsf_colon_rec$yvar$time, rsf_colon_rec$yvar$status), ignore_attr = TRUE)
    expect_equal(rsf_explainer$data, colon[, c("rx", "sex", "age", "obstruct", "perfor", "adhere", "nodes", "differ", "extent", "surg", "node4")], ignore_attr = TRUE)
})



test_that("default methods for creating explainers work correctly", {
    veteran <- survival::veteran

    ### survival::coxph ###

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)

    cph_exp <- explain(cph, verbose = FALSE)
    expect_s3_class(cph_exp, c("surv_explainer", "explainer"))
    expect_equal(cph_exp$label, "coxph")

    cph_without_params_1 <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran)
    expect_error(explain(cph_wihtout_params_1, verbose = FALSE))

    cph_without_params_2 <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE)
    expect_error(explain(cph_wihtout_params_2, verbose = FALSE))

    expect_output(explain(cph))


    ### ranger::ranger ###

    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)

    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    expect_s3_class(rsf_ranger_exp, c("surv_explainer", "explainer"))
    expect_equal(rsf_ranger_exp$label, "ranger")


    ### randomForestSRC::rfsrc ###

    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    rsf_src_exp <- explain(rsf_src, verbose = FALSE)
    expect_s3_class(rsf_src_exp, c("surv_explainer", "explainer"))
    expect_equal(rsf_src_exp$label, "rfsrc", ignore_attr = TRUE)


    ### rms::cph ###

    library(rms, quietly = TRUE)
    surv <- survival::Surv(veteran$time, veteran$status)
    cph <- rms::cph(surv ~ trt + celltype + karno + diagtime + age + prior,
                    data = veteran, surv=TRUE, model=TRUE, x=TRUE, y=TRUE)
    cph_rms_exp <- explain(cph, verbose = FALSE)
    expect_s3_class(cph_rms_exp, c("surv_explainer", "explainer"))
    expect_equal(cph_rms_exp$label, "coxph", ignore_attr = TRUE)


    ### parsnip::boost_tree ###

    library(censored, quietly = TRUE)
    bt <- parsnip::boost_tree() %>%
          parsnip::set_engine("mboost") %>%
          parsnip::set_mode("censored regression") %>%
          generics::fit(survival::Surv(time, status) ~ ., data = veteran)

    bt_exp <- explain(bt, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    expect_s3_class(bt_exp, c("surv_explainer", "explainer"))
    expect_equal(bt_exp$label, "model_fit_blackboost", ignore_attr = TRUE)
    detach("package:censored", unload = TRUE)

    ### explain.default ###

    exp_def <- explain.default(cph, verbose = FALSE)

    expect_s3_class(exp_def, c("surv_explainer", "explainer"))
    expect_equal(exp_def$label, "coxph")

    ## does DALEX::explain work ###

    exp_dal <- explain.default(rsf_src, verbose = FALSE)
    expect_s3_class(exp_dal, "explainer")
})


test_that("warnings in explain_survival work correctly", {

    veteran <- survival::veteran
    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    cph_exp <- explain(cph, verbose = FALSE, colorize = FALSE)

    expect_warning(explain_survival(cph,
                                    data = veteran,
                                    survival::Surv(veteran$time, veteran$status),
                                    predict_survival_function = pec::predictSurvProb,
                                    verbose = FALSE,
                                    label = matrix(1)))


    expect_warning(explain_survival(cph,
                                    data = veteran,
                                    times = c(1,2,3),
                                    predict_cumulative_hazard_function = pec::predictSurvProb,
                                    verbose = FALSE))
    expect_warning(
    expect_warning(explain_survival(cph,
                                    data = veteran,
                                    predict_survival_function = pec::predictSurvProb,
                                    verbose = FALSE,
                                    label = "my label")))

    expect_warning(explain_survival(cph,
                                    data = veteran,
                                    survival::Surv(veteran$time, veteran$status),
                                    times = c(1,2,3),
                                    predict_cumulative_hazard_function = pec::predictSurvProb,
                                    predict_function = "aaa",
                                    verbose = FALSE))

    expect_warning(explain_survival(cph,
                                    data = veteran,
                                    survival::Surv(veteran$time, veteran$status),
                                    times = c(1,2,3),
                                    predict_cumulative_hazard_function = pec::predictSurvProb,
                                    predict_survival_function = "aaa",
                                    verbose = FALSE))

    expect_warning(explain_survival(cph,
                                    data = veteran,
                                    survival::Surv(veteran$time, veteran$status),
                                    times = c(1,2,3),
                                    predict_survival_function = pec::predictSurvProb,
                                    predict_cumulative_hazard_function = "",
                                    verbose = FALSE))
})

test_that("default method for creating explainers for mlr3proba works correctly", {
    skip_on_ci()
    skip_on_cran()
    if (requireNamespace("mlr3proba", quietly=TRUE)){
        library(mlr3proba, quietly = TRUE)

        task <- TaskSurv$new("veteran", backend = veteran,
                             time = "time", event = "status")
        learner <- lrn("surv.coxph")
        learner$train(task)

        mlr3_exp <- explain(learner, y = task$truth(), data = as.data.frame(task$data())[, task$feature_names], verbose = FALSE)

        expect_s3_class(mlr3_exp, c("surv_explainer", "explainer"))
        expect_equal(mlr3_exp$label, "LearnerSurvCoxPH")

        detach("package:mlr3proba", unload = TRUE)
    }
})
