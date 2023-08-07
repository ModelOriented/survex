test_that("model_profile_2d with type = 'partial' works", {
    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_exp <- explain(rsf, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)

    mp_cph_pdp <- model_profile_2d(cph_exp,
                                   variable_splits_type = "uniform",
                                   variables = list(c("trt", "age"),
                                                    c("karno", "trt"),
                                                    c("karno", "age")),
                                   categorical_variables = "trt",
                                   grid_points = 6)
    mp_small <- model_profile_2d(cph_exp,
                                 variables = list(c("trt", "age")),
                                 categorical_variables = 1,
                                 grid_points = 2,
                                 N = 2)
    plot(mp_cph_pdp, times=cph_exp$times[1])

    mp_rsf_pdp <-  model_profile_2d(rsf_exp,
                                    variables = list(c("karno", "age")),
                                    grid_points = 6,
                                    output_type = "survival",
                                    N = 25)
    plot(mp_cph_pdp, mp_rsf_pdp, variables = list(c("karno", "age")), times=cph_exp$times[1])

    expect_output(print(mp_cph_pdp))
    expect_s3_class(mp_cph_pdp, "model_profile_2d_survival")
    expect_true(all(mp_cph_pdp$eval_times == cph_exp$times))
    expect_equal(ncol(mp_cph_pdp$result), 9)
    expect_true(all(unique(c(mp_cph_pdp$result$`_v1name_`, mp_cph_pdp$result$`_v2name_`))
                    %in% colnames(cph_exp$data)))

    expect_warning(plot(mp_cph_pdp))
    expect_error(model_profile_2d(rsf_exp))
    expect_error(model_profile_2d(rsf_exp, type = "conditional",
                                  variables = list(c("karno", "age"))))
    expect_error(model_profile_2d(rsf_exp, output_type = "risk",
                                  variables = list(c("karno", "age"))))
    }
)

test_that("model_profile_2d with type = 'accumulated' works", {
    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    rsf <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_exp <- explain(rsf, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)

    mp_rsf_ale <- model_profile_2d(rsf_exp,
                                variable_splits_type = "quantiles",
                                variables = list(c("karno", "age")),
                                grid_points = 6,
                                output_type = "survival",
                                type = "accumulated")

    mp_rsf_ale_noncentered <- model_profile_2d(rsf_exp,
                                   variable_splits_type = "quantiles",
                                   variables = list(c("karno", "age")),
                                   grid_points = 6,
                                   output_type = "survival",
                                   type = "accumulated",
                                   center = FALSE)

    expect_error(model_profile_2d(rsf_exp,
                                  type = "accumulated",
                                  variables = list(c("karno", "celltype"))))

    plot(mp_rsf_ale, times=rsf_exp$times[1])

    expect_output(print(mp_rsf_ale))
    expect_s3_class(mp_rsf_ale, "model_profile_2d_survival")
    expect_true(all(unique(mp_rsf_ale$eval_times) == rsf_exp$times))
    expect_equal(ncol(mp_rsf_ale$result), 14)
    expect_true(all(unique(c(mp_rsf_ale$result$`_v1name_`, mp_rsf_ale$result$`_v2name_`))
                    %in% colnames(rsf_exp$data)))

    expect_error(plot(mp_rsf_ale, variables = "nonexistent"))
    expect_error(plot(mp_rsf_ale,
                      variables = list(c("karno", "trt")),
                      categorical_variables="trt"))
    expect_error(model_profile_2d(mp_rsf_ale,
                      variables = list(c("karno", "trt")),
                      categorical_variables="trt",
                      type = "accumulated"))
    expect_error(plot(mp_rsf_ale, times = -1))
})
