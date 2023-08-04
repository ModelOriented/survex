test_that("model_performance works", {

    rotterdam <- survival::rotterdam
    rotterdam$pid <- NULL
    rotterdam <- rotterdam[sample(1:nrow(rotterdam), 300), ]

    cox_rotterdam_rec <- survival::coxph(survival::Surv(rtime, recur) ~ ., data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")], model = TRUE, x = TRUE, y = TRUE)
    rsf_rotterdam_rec <- ranger::ranger(survival::Surv(rtime, recur) ~ . - year - dtime - death - size, data = rotterdam, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 4)

    cph_exp_rot <- explain(cox_rotterdam_rec, y = survival::Surv(rotterdam$rtime, rotterdam$recur),  verbose = FALSE)
    rsf_exp_rot <- explain(rsf_rotterdam_rec, y = survival::Surv(rotterdam$rtime, rotterdam$recur), data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death", "rtime", "recur")], verbose = FALSE)

    cph_rot_perf <- model_performance(cph_exp_rot)
    rsf_rot_perf <- model_performance(rsf_exp_rot)

    expect_s3_class(cph_rot_perf, "list")
    expect_s3_class(rsf_rot_perf, "list")

    expect_lte(cph_rot_perf$result$`C-index`, 1)
    expect_gte(cph_rot_perf$result$`C-index`, 0)
    expect_lte(cph_rot_perf$result$`Integrated C/D AUC`, 1)
    expect_gte(cph_rot_perf$result$`Integrated C/D AUC`, 0)
    expect_lte(cph_rot_perf$result$`Integrated Brier score`, 1)
    expect_gte(cph_rot_perf$result$`Integrated Brier score`, 0)
    expect_equal(length(cph_rot_perf$result$`Brier score`), length(cph_exp_rot$times))
    expect_equal(length(cph_rot_perf$result$`C/D AUC`), length(cph_exp_rot$times))
    expect_true(all(cph_rot_perf$eval_times == cph_exp_rot$times))

    expect_lte(rsf_rot_perf$result$`C-index`, 1)
    expect_gte(rsf_rot_perf$result$`C-index`, 0)
    expect_lte(rsf_rot_perf$result$`Integrated C/D AUC`, 1)
    expect_gte(rsf_rot_perf$result$`Integrated C/D AUC`, 0)
    expect_lte(rsf_rot_perf$result$`Integrated Brier score`, 1)
    expect_gte(rsf_rot_perf$result$`Integrated Brier score`, 0)
    expect_equal(length(rsf_rot_perf$result$`Brier score`), length(rsf_exp_rot$times))
    expect_equal(length(rsf_rot_perf$result$`C/D AUC`), length(rsf_exp_rot$times))
    expect_true(all(rsf_rot_perf$eval_times == rsf_exp_rot$times))

    expect_output(print(cph_rot_perf))
    expect_output(print(rsf_rot_perf))

    plot(cph_rot_perf, rsf_rot_perf)
    plot(cph_rot_perf, rsf_rot_perf, rug = "events")
    plot(cph_rot_perf, rsf_rot_perf, rug = "censors")
    plot(cph_rot_perf, rsf_rot_perf, rug = "none")
    plot(cph_rot_perf, rsf_rot_perf, metrics_type = "scalar")

    cph_rot_perf_roc <- model_performance(cph_exp_rot, type = "roc", times = c(100, 200))
    rsf_rot_perf_roc <- model_performance(rsf_exp_rot, type = "roc", times = c(100, 200))

    plot(cph_rot_perf_roc)
    plot(rsf_rot_perf_roc)

    expect_error(model_performance(rsf_exp_rot, type = "roc", times = NULL))

})
