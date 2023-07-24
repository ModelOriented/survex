
veteran <- survival::veteran
rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran, respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)

test_that("global survshap explanations with kernelshap work for ranger, using new data", {

    ranger_global_survshap <- model_survshap(
        explainer = rsf_ranger_exp,
        new_observation = veteran[1:40, !colnames(veteran) %in% c("time", "status")],
        y_true = Surv(veteran$time[1:40], veteran$status[1:40]),
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap"
    )
    plot(ranger_global_survshap)

    expect_s3_class(ranger_global_survshap, c("aggregated_surv_shap", "surv_shap"))
    expect_equal(length(ranger_global_survshap$eval_times), length(rsf_ranger_exp$times))
    expect_true(all(names(ranger_global_survshap$variable_values) == colnames(rsf_ranger_exp$data)))

})

test_that("global survshap explanations with kernelshap work for ranger, using explainer data", {

    # using all explainer data
    ranger_global_survshap <- model_survshap(
        explainer = rsf_ranger_exp,
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap"
    )
    plot(ranger_global_survshap)

    # using only 6 observations
    ranger_global_survshap <- model_survshap(
        explainer = rsf_ranger_exp,
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap",
        N = 6L
    )
    plot(ranger_global_survshap)

})


test_that("global survshap explanations with treeshap work for ranger", {

    new_obs <- model.matrix(~ -1 + ., veteran[1:40, setdiff(colnames(veteran), c("time", "status"))])
    ranger_global_survshap_tree <- model_survshap(
        rsf_ranger_exp,
        new_obs,
        y_true = survival::Surv(veteran$time[1:40], veteran$status[1:40]),
        aggregation_method = "mean_absolute",
        calculation_method = "treeshap"
    )
    plot(ranger_global_survshap_tree)

    expect_s3_class(ranger_global_survshap_tree, c("aggregated_surv_shap", "surv_shap"))
    expect_equal(length(ranger_global_survshap_tree$eval_times), length(rsf_ranger_exp$times))
    expect_true(all(names(ranger_global_survshap_tree$variable_values) == colnames(rsf_ranger_exp$data)))

})
