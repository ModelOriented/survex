
# create objects here so that they do not have to be created redundantly
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
        new_observation = veteran[c(1:4, 17:20, 110:113, 126:129), !colnames(veteran) %in% c("time", "status")],
        y_true = Surv(veteran$time[c(1:4, 17:20, 110:113, 126:129)], veteran$status[c(1:4, 17:20, 110:113, 126:129)]),
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap"
    )
    ranger_global_survshap <- model_survshap(
        explainer = rsf_ranger_exp,
        new_observation = veteran[c(1:4, 17:20, 110:113, 126:129), !colnames(veteran) %in% c("time", "status")],
        y_true = Surv(veteran$time[c(1:4, 17:20, 110:113, 126:129)], veteran$status[c(1:4, 17:20, 110:113, 126:129)]),
        aggregation_method = "mean_absolute",
        calculation_method = "kernelshap",
        output_type = "chf"
    )
    plot(ranger_global_survshap)
    plot(ranger_global_survshap, geom = "beeswarm")
    plot(ranger_global_survshap, geom = "profile", variable = "karno", color_variable = "celltype")
    plot(ranger_global_survshap, geom = "profile", variable = "karno", color_variable = "age")
    plot(ranger_global_survshap, geom = "curves", variable = "karno")
    plot(ranger_global_survshap, geom = "curves", variable = "celltype")
    expect_output(plot(ranger_global_survshap, geom = "curves", variable = "karno", boxplot = TRUE))
    expect_error(plot(ranger_global_survshap, geom = "nonexistent"))

    single_survshap <- extract_predict_survshap(ranger_global_survshap, 5)
    expect_s3_class(single_survshap, c("predict_parts_survival", "surv_shap"))
    expect_error(extract_predict_survshap(ranger_global_survshap, 200))
    expect_error(extract_predict_survshap(single_survshap, 5))


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
    plot(cph_global_survshap, geom = "beeswarm")
    plot(cph_global_survshap, geom = "profile", variable = "karno", color_variable = "celltype")
    plot(cph_global_survshap, geom = "profile", variable = "karno", color_variable = "age")
    plot(cph_global_survshap, geom = "curves", variable = "karno")
    plot(cph_global_survshap, geom = "curves", variable = "celltype")
    expect_output(plot(cph_global_survshap, geom = "curves", variable = "karno", boxplot = TRUE))

    expect_s3_class(cph_global_survshap, c("aggregated_surv_shap", "surv_shap"))
    expect_equal(length(cph_global_survshap$eval_times), length(cph_exp$times))
    expect_true(all(names(cph_global_survshap$variable_values) == colnames(cph_exp$data)))
})

# testing if matrix works as input
rsf_ranger_matrix <- ranger::ranger(survival::Surv(time, status) ~ ., data = model.matrix(~ -1 + ., veteran), respect.unordered.factors = TRUE, num.trees = 100, mtry = 3, max.depth = 5)
rsf_ranger_exp_matrix <- explain(rsf_ranger_matrix, data = model.matrix(~ -1 + ., veteran[, -c(3, 4)]), y = survival::Surv(veteran$time, veteran$status), verbose = FALSE)

test_that("global survshap explanations with treeshap work for ranger", {

    new_obs <- model.matrix(~ -1 + ., veteran[1:40, setdiff(colnames(veteran), c("time", "status"))])
    ranger_global_survshap_tree <- model_survshap(
        rsf_ranger_exp_matrix,
        new_observation = new_obs,
        y_true = survival::Surv(veteran$time[1:40], veteran$status[1:40]),
        aggregation_method = "mean_absolute",
        calculation_method = "treeshap"
    )
    plot(ranger_global_survshap_tree)

    expect_s3_class(ranger_global_survshap_tree, c("aggregated_surv_shap", "surv_shap"))
    expect_equal(length(ranger_global_survshap_tree$eval_times), length(rsf_ranger_exp$times))
    expect_true(all(names(ranger_global_survshap_tree$variable_values) == colnames(rsf_ranger_exp_matrix$data)))

})
