test_that("survshap explanations work", {
    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = survival::Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    parts_cph <- predict_parts(cph_exp, veteran[1, !colnames(veteran) %in% c("time", "status")], y_true = matrix(c(100, 1), ncol = 2), aggregation_method = "sum_of_squares")
    parts_cph <- predict_parts(cph_exp, veteran[1, !colnames(veteran) %in% c("time", "status")], y_true = matrix(c(100, 1), ncol = 2), calculation_method = "exact_kernel", aggregation_method = "max_absolute")
    plot(parts_cph)
    plot(parts_cph, rug = "events")
    plot(parts_cph, rug = "censors")
    plot(parts_cph, rug = "none")

    parts_ranger <- predict_parts(rsf_ranger_exp, veteran[2, !colnames(veteran) %in% c("time", "status")], y_true = c(100, 1), aggregation_method = "mean_absolute")
    plot(parts_ranger)

    # test ranger with kernelshap when using a matrix as input for data and new observation
    rsf_ranger_matrix <- ranger::ranger(survival::Surv(time, status) ~ ., data = model.matrix(~ -1 + ., veteran), respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_ranger_exp_matrix <- explain(rsf_ranger_matrix, data = model.matrix(~ -1 + ., veteran[, -c(3, 4)]), y = survival::Surv(veteran$time, veteran$status), verbose = FALSE)
    new_obs <- model.matrix(~ -1 + ., veteran[2, !colnames(veteran) %in% c("time", "status")])
    parts_ranger_kernelshap <- predict_parts(
        rsf_ranger_exp_matrix,
        new_observation = new_obs,
        y_true = c(100, 1),
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap"
    )
    plot(parts_ranger_kernelshap)

    parts_src <- predict_parts(rsf_src_exp, veteran[3, !colnames(veteran) %in% c("time", "status")])
    plot(parts_src)

    plot(parts_cph, parts_ranger, parts_src)

    parts_cph2 <- predict_parts(cph_exp, veteran[4,], explanation_label = "second_explanation")
    plot(parts_cph, parts_cph2)

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

test_that("local survshap explanations with treeshap work for ranger", {

    veteran <- survival::veteran

    rsf_ranger_matrix <- ranger::ranger(survival::Surv(time, status) ~ ., data = model.matrix(~ -1 + ., veteran), respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_ranger_exp_matrix <- explain(rsf_ranger_matrix, data = model.matrix(~ -1 + ., veteran[, -c(3, 4)]), y = survival::Surv(veteran$time, veteran$status), verbose = FALSE)


    new_obs <- model.matrix(~ -1 + ., veteran[2, setdiff(colnames(veteran), c("time", "status"))])
    parts_ranger <- model_survshap(
        rsf_ranger_exp_matrix,
        new_obs,
        y_true = c(veteran$time[2], veteran$status[2]),
        aggregation_method = "mean_absolute",
        calculation_method = "treeshap"
    )
    plot(parts_ranger)

    expect_s3_class(parts_ranger, c("predict_parts_survival", "surv_shap"))
    expect_equal(nrow(parts_ranger$result), length(rsf_ranger_exp_matrix$times))
    expect_true(all(colnames(parts_ranger$result) == colnames(rsf_ranger_exp_matrix$data)))

})

test_that("survshap explanations with output_type = 'chf' work", {
    veteran <- survival::veteran

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    parts_cph <- predict_parts(cph_exp, veteran[1, !colnames(veteran) %in% c("time", "status")], y_true = matrix(c(100, 1), ncol = 2), aggregation_method = "sum_of_squares", output_type = "chf")
    parts_cph <- predict_parts(cph_exp, veteran[1, !colnames(veteran) %in% c("time", "status")], y_true = matrix(c(100, 1), ncol = 2), calculation_method = "exact_kernel", aggregation_method = "max_absolute", output_type = "chf")
    plot(parts_cph)
    plot(parts_cph, rug = "events")
    plot(parts_cph, rug = "censors")
    plot(parts_cph, rug = "none")

    parts_ranger <- predict_parts(rsf_ranger_exp, veteran[2, !colnames(veteran) %in% c("time", "status")], y_true = c(100, 1), aggregation_method = "mean_absolute", output_type = "chf")
    plot(parts_ranger)

    parts_src <- predict_parts(rsf_src_exp, veteran[3, !colnames(veteran) %in% c("time", "status")], output_type = "chf")
    plot(parts_src)

    plot(parts_cph, parts_ranger, parts_src)

    parts_cph2 <- predict_parts(cph_exp, veteran[4,], explanation_label = "second_explanation", output_type = "chf")
    plot(parts_cph, parts_cph2)

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
    ranger_survlime <- predict_parts(rsf_ranger_exp, new_observation = veteran[1, -c(3, 4)], type = "survlime", sample_around_instance = FALSE)
    rsf_survlime <- predict_parts(rsf_src_exp, new_observation = veteran[1, -c(3, 4)], type = "survlime", categorical_variables = 1)

    # error on to few columns
    expect_error(predict_parts(rsf_src_exp, new_observation = veteran[1, -c(1, 2 ,3, 4)], type = "survlime"))

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
