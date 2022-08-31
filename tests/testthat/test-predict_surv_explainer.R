test_that("prediction funciton works on explainer", {
    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    cox_risk <- predict(cph_exp, veteran[1:10, ], output_type = "risk")
    cox_sf   <- predict(cph_exp, veteran[1:10, ], output_type = "survival")
    cox_chf  <- predict(cph_exp, veteran[1:10, ], output_type = "chf")

    expect_length(cox_risk, 10)
    expect_equal(ncol(cox_sf), length(cph_exp$times))
    expect_equal(nrow(cox_sf), 10)
    expect_equal(ncol(cox_chf), length(cph_exp$times))
    expect_equal(nrow(cox_chf), 10)
    expect_true(all(cox_sf <= 1))
    expect_true(all(cox_sf >= 0))
    expect_true(all(cox_risk >= 0))
    expect_true(all(cox_chf >= 0))

    ranger_risk <- predict(rsf_ranger_exp, veteran[1:10, ], output_type = "risk")
    ranger_sf   <- predict(rsf_ranger_exp, veteran[1:10, ], output_type = "survival")
    ranger_chf  <- predict(rsf_ranger_exp, veteran[1:10, ], output_type = "chf")

    expect_length(ranger_risk, 10)
    expect_equal(ncol(ranger_sf), length(rsf_ranger_exp$times))
    expect_equal(nrow(ranger_sf), 10)
    expect_equal(ncol(ranger_chf), length(rsf_ranger_exp$times))
    expect_equal(nrow(ranger_chf), 10)
    expect_true(all(ranger_sf <= 1))
    expect_true(all(ranger_sf >= 0))
    expect_true(all(ranger_risk >= 0))
    expect_true(all(ranger_chf >= 0))

    src_risk <- predict(rsf_src_exp, veteran[1:10, ], output_type = "risk")
    src_sf   <- predict(rsf_src_exp, veteran[1:10, ], output_type = "survival")
    src_chf  <- predict(rsf_src_exp, veteran[1:10, ], output_type = "chf")

    expect_length(src_risk, 10)
    expect_equal(ncol(src_sf), length(rsf_src_exp$times))
    expect_equal(nrow(src_sf), 10)
    expect_equal(ncol(src_chf), length(rsf_src_exp$times))
    expect_equal(nrow(src_chf), 10)
    expect_true(all(src_sf <= 1))
    expect_true(all(src_sf >= 0))
    expect_true(all(src_risk >= 0))
    expect_true(all(src_chf >= 0))

    expect_output(print(predict(cph_exp)))
    expect_error(predict(cph_exp, output_type = "non_existent"))

})
