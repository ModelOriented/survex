test_that("C-index fpi works", {

    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL

    rotterdam <- rotterdam[c(1:15, 2001:2015), ]

    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )

    rsf_rotterdam_rec <- ranger::ranger(survival::Surv(rotterdam$rtime, rotterdam$recur) ~ .,
                                        data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
                                        respect.unordered.factors = TRUE,
                                        num.trees = 100,
                                        mtry = 3,
                                        max.depth = 4)


    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]
    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    cph_exp <- explain(cph, verbose = FALSE)


    coxph_explainer <- explain(cox_rotterdam_rec, y = survival::Surv(rotterdam$rtime, rotterdam$recur), verbose = FALSE)
    forest_explainer <- explain(rsf_rotterdam_rec,
                                data = rotterdam[,  !colnames(rotterdam) %in% c("year", "dtime", "death", "rtime", "recur")],
                                y = survival::Surv(rotterdam$rtime, rotterdam$recur), verbose = FALSE)


    mp_cph_cind <- model_parts(coxph_explainer, loss = loss_one_minus_c_index, type = "variable_importance", output_type = "risk")
    mp_rsf_cind <- model_parts(forest_explainer, loss = loss_one_minus_c_index, output_type = "risk")

    expect_equal(nrow(mp_rsf_cind[mp_rsf_cind$permutation == 0, ]), ncol(forest_explainer$data) + 2)
    expect_equal(nrow(mp_cph_cind[mp_cph_cind$permutation == 0, ]), ncol(coxph_explainer$data) + 2)
    expect_s3_class(mp_cph_cind, "model_parts")


    cph_model_parts_dalex <- model_parts(cph_exp, loss_function = loss_one_minus_c_index,
                                         output_type = "risk")

    expect_true(all(cph_model_parts_dalex$dropout_loss <= 1))
    expect_true(all(cph_model_parts_dalex$dropout_loss > 0))
    expect_s3_class(cph_model_parts_dalex, "model_parts")

    expect_error(model_parts(coxph_explainer, type = "nonexistent"))
    expect_error(model_parts(coxph_explainer, output_type = "nonexistent"))

    plot(cph_model_parts_dalex)

})



test_that("Brier score fpi works", {

    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    cph_model_parts_brier <- model_parts(cph_exp)
    rsf_ranger_model_parts <- model_parts(rsf_ranger_exp)

    expect_s3_class(cph_model_parts_brier, "model_parts_survival")
    expect_s3_class(rsf_ranger_model_parts, "model_parts_survival")
    expect_equal(ncol(cph_model_parts_brier$result), ncol(cph_exp$data) + 5) # times, full_model, permutation, baseline, label
    expect_equal(nrow(cph_model_parts_brier$result[cph_model_parts_brier$result$`_permutation_` == 0, ]), length(cph_exp$times))

    plot(cph_model_parts_brier)
    plot(cph_model_parts_brier, rsf_ranger_model_parts, max_vars = 3)

    # specifying loss function brier
    rsf_src_model_parts_brier <- model_parts(rsf_src_exp, loss_function = loss_brier_score, output_type = "survival")
    plot(rsf_src_model_parts_brier)

    expect_s3_class(rsf_src_model_parts_brier, "model_parts_survival")
    expect_equal(ncol(rsf_src_model_parts_brier$result), ncol(cph_exp$data) + 5) # times, full_model, permutation, baseline, label

    expect_output(print(cph_model_parts_brier))


    ### groups
    mp_groups_1 <- model_parts(cph_exp, type = "ratio", variable_groups = list(group1 = c("celltype", "trt"), group2 = c("age", "prior")))
    mp_groups_2 <- model_parts(cph_exp, type = "difference", variable_groups = list(group1 = c("celltype", "trt"), group2 = c("age", "prior")), N = 70)
    plot(mp_groups_2)

    expect_error(model_parts(cph_exp, variable_groups = list(group1 = c("sss", "ss"), group2 = c("f", "f"))))
    expect_error(model_parts(cph_exp, variable_groups = list(group1 = c("sss", "ss"), group2 = c("f", "f"))))
    expect_warning(model_parts(cph_exp, variable_groups = list(c("celltype", "trt"), c("age", "prior"))))

})


test_that("CD/AUC fpi works", {

    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    cph_model_parts_auc <- model_parts(cph_exp, loss = loss_one_minus_cd_auc, B = 1)
    rsf_ranger_model_auc <- model_parts(rsf_ranger_exp, loss = loss_one_minus_cd_auc, B = 1)

    expect_s3_class(cph_model_parts_auc, "model_parts_survival")
    expect_s3_class(rsf_ranger_model_auc, "model_parts_survival")
    expect_equal(ncol(cph_model_parts_auc$result), ncol(cph_exp$data) + 5) # times, full_model, permutation, baseline, label

    plot(cph_model_parts_auc)
    plot(cph_model_parts_auc, rsf_ranger_model_auc)

})



test_that("integrated metrics fpi works", {

    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
    rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = veteran)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
    rsf_src_exp <- explain(rsf_src, verbose = FALSE)

    # auc
    cph_model_parts_int_auc <- model_parts(cph_exp, loss = loss_one_minus_integrated_cd_auc, B = 1)
    rsf_ranger_model_parts_int_auc <- model_parts(rsf_ranger_exp, loss = loss_one_minus_integrated_cd_auc, B = 1)

    expect_equal(nrow(cph_model_parts_int_auc[cph_model_parts_int_auc$permutation == 0, ]), ncol(cph_exp$data) + 2)
    expect_s3_class(cph_model_parts_int_auc, "model_parts_survival")
    expect_true(all(cph_model_parts_int_auc$dropout_loss <= 1))
    expect_true(all(cph_model_parts_int_auc$dropout_loss >= 0))

    expect_equal(nrow(rsf_ranger_model_parts_int_auc[rsf_ranger_model_parts_int_auc$permutation == 0, ]), ncol(rsf_ranger_exp$data) + 2)
    expect_s3_class(rsf_ranger_model_parts_int_auc, "model_parts_survival")
    expect_true(all(rsf_ranger_model_parts_int_auc$dropout_loss <= 1))
    expect_true(all(rsf_ranger_model_parts_int_auc$dropout_loss >= 0))

    plot(cph_model_parts_int_auc, rsf_ranger_model_parts_int_auc)

    # brier

    cph_model_parts_int_brier <- model_parts(cph_exp, loss = loss_integrated_brier_score, B = 10)
    rsf_ranger_model_parts_int_brier <- model_parts(rsf_ranger_exp, loss = loss_integrated_brier_score, B = 10)

    expect_equal(nrow(cph_model_parts_int_brier[cph_model_parts_int_brier$permutation == 0, ]), ncol(cph_exp$data) + 2)
    expect_s3_class(cph_model_parts_int_brier, "model_parts_survival")
    expect_true(all(cph_model_parts_int_brier$dropout_loss <= 1))
    expect_true(all(cph_model_parts_int_brier$dropout_loss >= 0))

    expect_equal(nrow(rsf_ranger_model_parts_int_brier[rsf_ranger_model_parts_int_brier$permutation == 0, ]), ncol(rsf_ranger_exp$data) + 2)
    expect_s3_class(rsf_ranger_model_parts_int_brier, "model_parts_survival")
    expect_true(all(rsf_ranger_model_parts_int_brier$dropout_loss <= 1))
    expect_true(all(rsf_ranger_model_parts_int_brier$dropout_loss >= 0))

    plot(cph_model_parts_int_brier, rsf_ranger_model_parts_int_brier)


    expect_error(model_parts(list(data = NULL), loss = loss_integrated_brier_score))
    expect_error(model_parts(list(data = "a", y = NULL), loss = loss_integrated_brier_score))
    expect_error(model_parts(list(data = NULL), loss = loss_integrated_brier_score, variable_groups = "a"))

    expect_error(model_parts(cph_exp, loss = loss_integrated_brier_score, variable_groups = list(group1 = c("sss", "ss"), group2 = c("f", "f"))))
    expect_error(model_parts(cph_exp, loss = loss_integrated_brier_score, variable_groups = list(group1 = c("sss", "ss"), group2 = c("f", "f"))))
    expect_warning(model_parts(cph_exp, loss = loss_integrated_brier_score, variable_groups = list(c("celltype", "trt"), c("age", "prior"))))

    mp_groups_1 <- model_parts(cph_exp, loss = loss_integrated_brier_score, type = "ratio", variable_groups = list(group1 = c("celltype", "trt"), group2 = c("age", "prior")))
    mp_groups_2 <- model_parts(cph_exp, loss = loss_integrated_brier_score, type = "difference", variable_groups = list(group1 = c("celltype", "trt"), group2 = c("age", "prior")), N = 70)
    plot(mp_groups_2)

})
