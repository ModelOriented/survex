test_that("survshap explanations work", {
    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    parts_cph <- predict_parts(cph_exp, veteran[1, !colnames(veteran) %in% c("time", "status")], y_true = matrix(c(100, 1), ncol = 2), aggregation_method = "sum_of_squares")
    plot(parts_cph)

    parts_ranger <- predict_parts(rsf_ranger_exp, veteran[2, !colnames(veteran) %in% c("time", "status")], y_true = c(100, 1), aggregation_method = "mean_absolute")
    plot(parts_ranger)

    parts_src <- predict_parts(rsf_src_exp, veteran[3, !colnames(veteran) %in% c("time", "status")])
    plot(parts_src)

    plot(parts_cph, parts_ranger, parts_src)

    expect_s3_class(parts_cph, c("predict_parts_survival", "surv_shap"))
    expect_s3_class(parts_ranger, c("predict_parts_survival", "surv_shap"))
    expect_s3_class(parts_src, c("predict_parts_survival", "surv_shap"))

    expect_equal(nrow(parts_cph$result), length(cph_exp$times))
    expect_equal(nrow(parts_ranger$result), length(rsf_ranger_exp$times))
    expect_equal(nrow(parts_src$result), length(rsf_src_exp$times))

    expect_true(all(colnames(parts_cph$result) == colnames(cph_exp$data)))
    expect_true(all(colnames(parts_ranger$result) == colnames(rsf_ranger_exp$data)))
    expect_true(all(colnames(parts_src$result) == colnames(rsf_src_exp$data)))

    expect_output(print(parts_cph))

    expect_error(predict_parts(cph_exp, veteran[1, ], aggregation_method = "nonexistent"))
    expect_error(predict_parts(cph_exp, veteran[1, ], calculation_method = "sampling"))
    expect_error(predict_parts(cph_exp, veteran[1, ], calculation_method = "nonexistent"))
    expect_error(predict_parts(cph_exp, veteran[1, c(1, 1, 1, 1, 1)], calculation_method = "nonexistent"))


})


test_that("survlime explanations work", {

    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    cph_survlime <- predict_parts(cph_exp, new_observation = veteran[1, -c(3, 4)], type = "survlime")
    ranger_survlime <- predict_parts(rsf_ranger_exp, new_observation = veteran[1, -c(3, 4)], type = "survlime")
    rsf_survlime <- predict_parts(rsf_src_exp, new_observation = veteran[1, -c(3, 4)], type = "survlime")

    plot(cph_survlime, type = "coefficients")
    plot(cph_survlime, type = "local_importance")
    plot(cph_survlime, show_survival_function = FALSE)

    expect_error(plot(cph_survlime, type = "nonexistent"))

    expect_s3_class(cph_survlime, c("predict_parts_survival", "surv_lime"))
    expect_s3_class(ranger_survlime, c("predict_parts_survival", "surv_lime"))
    expect_s3_class(rsf_survlime, c("predict_parts_survival", "surv_lime"))

    expect_gte(length(cph_survlime$result), ncol(cph_exp$data))
    expect_gte(length(ranger_survlime$result), ncol(rsf_ranger_exp$data))
    expect_gte(length(rsf_survlime$result), ncol(rsf_src_exp$data))

    expect_setequal(cph_survlime$black_box_sf_times, cph_exp$times)
    expect_setequal(ranger_survlime$black_box_sf_times, rsf_ranger_exp$times)
    expect_setequal(rsf_survlime$black_box_sf_times, rsf_src_exp$times)

    expect_output(print(cph_survlime))

})

test_that("default DALEX::predict_parts is ok", {

    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)

    cph_exp <- explain(cph, verbose = FALSE)

    cph_pp <- predict_parts(cph_exp, veteran[2, -c(3, 4)], output_type = "risk", type = "shap")

    expect_output(print(cph_pp))
})
