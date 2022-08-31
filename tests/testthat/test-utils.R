test_that("changing survival to cumulative hazard works", {

    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL
    rsf_rotterdam_rec <- ranger::ranger(survival::Surv(rtime, recur) ~ . - year - dtime - death,
                                        data = rotterdam,
                                        respect.unordered.factors = TRUE,
                                        num.trees = 100,
                                        mtry = 3,
                                        max.depth = 4)


    colon <- survival::colon
    colon$study <- NULL
    colon$id <- NULL
    colon$extent <- as.factor(colon$extent)
    colon <- na.omit(colon)
    rsf_colon_rec <- randomForestSRC::rfsrc(Surv(time, status) ~ rx + sex + age + obstruct + perfor +
                                                adhere + nodes + differ + extent + surg + node4,
                                            data = colon, subset = colon$etype == 1)


    # single prediction
    prediction1_ranger <- predict(rsf_rotterdam_rec, rotterdam[1, ])
    prediction1_rfsrc <- predict(rsf_colon_rec, colon[42, ])

    expect_equal(survival_to_cumulative_hazard(prediction1_ranger$survival[1:10]), prediction1_ranger$chf[1:10])
    expect_equal(survival_to_cumulative_hazard(prediction1_rfsrc$survival[1:10]), prediction1_rfsrc$chf[1:10], tolerance = 0.05)
    # multiple predictions
    prediction2_ranger <- predict(rsf_rotterdam_rec, rotterdam[c(1, 2, 3, 4), ])
    prediction2_rfsrc <- predict(rsf_colon_rec, colon[c(1, 2, 3, 4), ])
    expect_equal(survival_to_cumulative_hazard(prediction2_ranger$survival[, 1:100]), prediction2_ranger$chf[, 1:100])
    expect_equal(survival_to_cumulative_hazard(prediction2_rfsrc$survival[, 1:100]), prediction2_rfsrc$chf[, 1:100], tolerance = 0.05)

})


test_that("changing cumulative hazard to survival works", {

    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL
    rsf_rotterdam_rec <- ranger::ranger(survival::Surv(rtime, recur) ~ . - year - dtime - death,
                                        data = rotterdam,
                                        respect.unordered.factors = TRUE,
                                        num.trees = 100,
                                        mtry = 3,
                                        max.depth = 4)
    colon <- survival::colon
    colon$study <- NULL
    colon$id <- NULL
    colon$extent <- as.factor(colon$extent)
    colon <- na.omit(colon)
    rsf_colon_rec <- randomForestSRC::rfsrc(Surv(time, status) ~ rx + sex + age + obstruct + perfor +
                                                adhere + nodes + differ + extent + surg + node4,
                                            data = colon, subset = colon$etype == 1)

    # single prediction
    prediction1_ranger <- predict(rsf_rotterdam_rec, rotterdam[1, ])
    prediction1_rfsrc <- predict(rsf_colon_rec, colon[42, ])

    expect_equal(cumulative_hazard_to_survival(prediction1_ranger$chf), prediction1_ranger$survival)
    expect_equal(cumulative_hazard_to_survival(prediction1_rfsrc$chf), prediction1_rfsrc$survival, tolerance = 0.05)
    # multiple predictions
    prediction2_ranger <- predict(rsf_rotterdam_rec, rotterdam[c(1, 2, 3, 4), ])
    prediction2_rfsrc <- predict(rsf_colon_rec, colon[c(1, 2, 3, 4), ])
    expect_equal(cumulative_hazard_to_survival(prediction2_ranger$chf), prediction2_ranger$survival)
    expect_equal(cumulative_hazard_to_survival(prediction2_rfsrc$chf), prediction2_rfsrc$survival, tolerance = 0.05)

})

test_that("explainer checkers work", {


    cph_exp <- suppressWarnings(explain_survival("cph", y=matrix(c(1,1,1,1,1,1), ncol=2), data=NULL, predict_function=NULL, verbose=FALSE))
    cph_exp$predict_function <- NULL
    cph_exp$y <- NULL

    expect_error(test_explainer("a", "test"))
    expect_error(test_explainer(cph_exp, "test", has_data = T))
    expect_error(test_explainer(cph_exp, "test", has_y = T))
    expect_error(test_explainer(cph_exp, "test", has_chf = T))
    expect_error(test_explainer(cph_exp, "test", has_survival = T))
    expect_error(test_explainer(cph_exp, "test", has_predict = T))

    expect_error(check_times("random words"))

})


