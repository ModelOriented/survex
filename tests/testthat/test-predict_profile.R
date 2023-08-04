test_that("ceteris_paribus works", {

    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    cph_pp <- predict_profile(cph_exp, veteran[2, -c(3, 4)])
    ranger_pp <- predict_profile(rsf_ranger_exp, veteran[2, -c(3, 4)])

    plot(cph_pp, colors = c("#ff0000", "#00ff00", "#0000ff"))
    plot(cph_pp)
    plot(cph_pp, numerical_plot_type = "contours")
    plot(cph_pp, ranger_pp, rug = "events", variables = c("karno", "age"))
    plot(cph_pp, rug = "censors", variable_type = "numerical")
    plot(cph_pp, rug = "none")

    expect_error(plot(cph_pp, variable_type = "nonexistent"))
    expect_error(plot(cph_pp, numerical_plot_type = "nonexistent"))

    cph_pp_cat <- predict_profile(cph_exp, veteran[2, -c(3, 4)], variables = c("celltype"))
    plot(cph_pp_cat, variable_type = "categorical", colors = c("#ff0000", "#00ff00", "#0000ff"))
    plot(cph_pp_cat, variable_type = "categorical")


    expect_s3_class(cph_pp, c("predict_profile_survival", "surv_ceteris_paribus"))
    expect_s3_class(cph_pp_cat, c("predict_profile_survival", "surv_ceteris_paribus"))

    expect_true(all(unique(cph_pp$result$`_vname_`) %in% colnames(cph_exp$data)))
    expect_true(all(unique(cph_pp_cat$result$`_vname_`) %in% colnames(cph_exp$data)))

    expect_true(all(unique(cph_pp$result$`_yhat_`) <= 1))
    expect_true(all(unique(cph_pp_cat$result$`_yhat_`) <= 1))

    expect_true(all(unique(cph_pp$result$`_yhat_`) >= 0))
    expect_true(all(unique(cph_pp_cat$result$`_yhat_`)  >= 0))

    expect_setequal(cph_pp$eval_times, cph_exp$times)
    expect_setequal(cph_pp_cat$eval_times, cph_exp$times)

    expect_output(print(cph_pp))

})


test_that("default DALEX::ceteris_paribus works", {

    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)

    cph_exp <- explain(cph, verbose = FALSE)

    cph_pp <- predict_profile(cph_exp, veteran[2, -c(3, 4)], output_type = "risk")

    expect_output(print(cph_pp))
})


test_that("error is thrown for incorrect type", {
    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)

    cph_exp <- explain(cph, verbose = FALSE)

    expect_error(predict_profile(cph_exp, veteran[2, -c(3, 4)], output_type = "survival", type = "nonexisting type"))
})
