test_that("model_profile with type = 'partial' works", {

    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)


    mp_cph_cat <- model_profile(cph_exp, output_type = "survival", variable_splits_type = "quantiles", variable_type = "categorical", grid_points = 6)
    plot(mp_cph_cat, variables = "celltype", variable_type = "categorical")

    expect_s3_class(mp_cph_cat, "model_profile_survival")
    expect_true(all(mp_cph_cat$eval_times == cph_exp$times))
    expect_equal(ncol(mp_cph_cat$result), 7)
    expect_true(all(unique(mp_cph_cat$result$`_vname_`) %in% colnames(cph_exp$data)))


    mp_cph_num <- model_profile(cph_exp, output_type = "survival", variable_splits_type = "quantiles", variable_type = "numerical", grid_points = 6)
    plot(mp_cph_num, variable_type = "numerical")
    plot(mp_cph_num, numerical_plot_type = "contours")

    expect_s3_class(mp_cph_num, "model_profile_survival")
    expect_true(all(unique(mp_cph_num$eval_times) == cph_exp$times))
    expect_equal(ncol(mp_cph_num$result), 7)
    expect_true(all(unique(mp_cph_num$result$`_vname_`) %in% colnames(cph_exp$data)))


    mp_rsf_cat <- model_profile(rsf_ranger_exp, output_type = "survival", variable_splits_type = "uniform", variable_type = "categorical", grid_points = 6)
    plot(mp_rsf_cat, variable_type = "categorical")

    expect_s3_class(mp_rsf_cat, "model_profile_survival")
    expect_true(all(mp_rsf_cat$eval_times == cph_exp$times))
    expect_equal(ncol(mp_rsf_cat$result), 7)
    expect_true(all(unique(mp_rsf_cat$result$`_vname_`) %in% colnames(rsf_ranger_exp$data)))


    mp_rsf_num <- model_profile(rsf_ranger_exp, output_type = "survival", variable_splits_type = "uniform", variable_type = "numerical", grid_points = 6)
    plot(mp_rsf_num, variable_type = "numerical")
    plot(mp_rsf_num, variable_type = "numerical", numerical_plot_type = "contours")

    expect_s3_class(mp_rsf_num, "model_profile_survival")
    expect_true(all(mp_rsf_num$eval_times == cph_exp$times))
    expect_equal(ncol(mp_rsf_num$result), 7)
    expect_true(all(unique(mp_rsf_num$result$`_vname_`) %in% colnames(rsf_ranger_exp$data)))

    expect_output(print(mp_cph_num))
    expect_error(plot(mp_rsf_num, variables = "nonexistent", grid_points = 6))


    model_profile(rsf_ranger_exp, varaibles = "trt", grid_points = 6)
    expect_error(model_profile(rsf_ranger_exp, type = "conditional"))

    })

test_that("model_profile with type = 'accumulated' works", {

    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]
    type <- 'accumulated'

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)


    mp_cph_cat <- model_profile(cph_exp,
                                output_type = "survival",
                                variable_type = "categorical",
                                grid_points = 6,
                                type = type)
    plot(mp_cph_cat, variables = "celltype", variable_type = "categorical")

    expect_s3_class(mp_cph_cat, "model_profile_survival")
    expect_true(all(mp_cph_cat$eval_times == cph_exp$times))
    expect_equal(ncol(mp_cph_cat$result), 7)
    expect_true(all(unique(mp_cph_cat$result$`_vname_`) %in% colnames(cph_exp$data)))


    mp_cph_num <- model_profile(cph_exp,
                                output_type = "survival",
                                variable_type = "numerical",
                                grid_points = 6,
                                type = type)
    plot(mp_cph_num, variable_type = "numerical")
    plot(mp_cph_num, numerical_plot_type = "contours")

    expect_s3_class(mp_cph_num, "model_profile_survival")
    expect_true(all(unique(mp_cph_num$eval_times) == cph_exp$times))
    expect_equal(ncol(mp_cph_num$result), 7)
    expect_true(all(unique(mp_cph_num$result$`_vname_`) %in% colnames(cph_exp$data)))

    expect_output(print(mp_cph_num))
    expect_error(plot(mp_rsf_num, variables = "nonexistent", grid_points = 6))
})

test_that("default DALEX::model_profile is ok", {

    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    cph_exp <- explain(cph, verbose = FALSE)
    cph_mp <- model_profile(cph_exp, output_type = "risk")

    expect_output(print(cph_mp))

    expect_error(model_profile(cph_exp, output_type = "something_else"))
})
