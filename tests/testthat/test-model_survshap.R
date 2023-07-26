
veteran <- survival::veteran
rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)

cph <- survival::coxph(
        survival::Surv(time, status) ~ .,
        data = veteran,
        model = TRUE,
        x = TRUE,
        y = TRUE
    )
cph_exp <- explain(cph, verbose = FALSE)

test_that("global survshap explanations with kernelshap work for ranger, using new data", {

    ranger_global_survshap <- model_survshap(
        explainer = rsf_ranger_exp,
        new_observation = veteran[1:10, !colnames(veteran) %in% c("time", "status")],
        y_true = Surv(veteran$time[1:10], veteran$status[1:10]),
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap"
    )
    plot(ranger_global_survshap)
    plot(ranger_global_survshap, kind = "swarm")
    plot(ranger_global_survshap, kind = "profile")
    plot(ranger_global_survshap, kind = "profile", variable = "karno", color_variable = "celltype")
    plot(ranger_global_survshap, kind = "profile", variable = "karno", color_variable = "age")


    expect_s3_class(ranger_global_survshap, c("aggregated_surv_shap", "surv_shap"))
    expect_equal(length(ranger_global_survshap$eval_times), length(rsf_ranger_exp$times))
    expect_true(all(names(ranger_global_survshap$variable_values) == colnames(rsf_ranger_exp$data)))
})

test_that("global survshap explanations with kernelshap work for coxph, using explainer data", {

    # using all explainer data
    cph_global_survshap <- model_survshap(
        explainer = cph_exp,
        calculation_method = "kernelshap"
    )
    plot(cph_global_survshap)
    plot(cph_global_survshap, kind = "swarm")
    plot(cph_global_survshap, kind = "profile")
    plot(cph_global_survshap, kind = "profile", variable = "karno", color_variable = "celltype")
    plot(cph_global_survshap, kind = "profile", variable = "karno", color_variable = "age")

    expect_s3_class(cph_global_survshap, c("aggregated_surv_shap", "surv_shap"))
    expect_equal(length(cph_global_survshap$eval_times), length(cph_exp$times))
    expect_true(all(names(cph_global_survshap$variable_values) == colnames(cph_exp$data)))
})
