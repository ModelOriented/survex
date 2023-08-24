test_that("model_diagnostics for survival residuals works", {
    veteran <- survival::veteran[c(1:3, 16:18, 46:48, 56:58, 71:73, 91:93, 111:113, 126:128), ]

    cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)
    rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)

    cph_exp <- explain(cph, verbose = FALSE)
    rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)

    md_cph <- model_diagnostics(cph_exp)
    md_rsf <- model_diagnostics(rsf_ranger_exp)
    expect_s3_class(md_cph, "model_diagnostics_survival")
    expect_true(all(md_cph$result$time == cph_exp$y[,1]))
    expect_equal(ncol(md_cph$result) - ncol(cph_exp$data), 6)

    plot(md_rsf)
    plot(md_rsf, plot_type = "martingale")
    plot(md_rsf, plot_type = "Cox-Snell")
    plot(md_cph, md_rsf, xvariable = "age")
    plot(md_cph, md_rsf, smooth = FALSE)

    expect_error(plot(md_cph, md_rsf, xvariable = "nonexistent"))
    expect_error(plot(md_cph, md_rsf, plot_type = "nonexistent"))
    expect_error(plot(md_cph, exp))
})
