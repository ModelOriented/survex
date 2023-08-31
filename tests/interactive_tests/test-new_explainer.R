test_that("all functionality of a new explainer works correctly", {

    ### DATASET ###
    veteran <- survival::veteran
    veteran_data <- veteran[, -c(3, 4)]
    veteran_y <- survival::Surv(veteran$time, veteran$status)

    ### MODEL CREATION ###

    library(mlr3proba)
    library(mlr3learners)
    library(mlr3extralearners)
    library(mlr3pipelines)

    veteran$celltype <- NULL
    task <- TaskSurv$new("veteran", backend = veteran,
                         time = "time", event = "status")

    #-- TESTS --#
    expect_s3_class(exp, c("surv_explainer", "explainer"))
    expect_true(!is.null(exp$label))

    #########################
    ### MODEL PERFORMANCE ###
    #########################

    m_perf <- model_performance(exp)
    m_perf_roc <- model_performance(exp, type = "roc", times = c(100, 300))



    plot(m_perf)
    plot(m_perf, metrics_type = "scalar")

    plot(m_perf_roc)

    #-- TESTS --#

    expect_s3_class(m_perf, c("surv_model_performance", "list"))
    expect_false(is.null(m_perf$eval_times))
    expect_false(is.null(m_perf$brier_score))
    expect_false(is.null(m_perf$auc))
    expect_false(is.null(m_perf$cindex))
    expect_false(is.null(m_perf$iauc))
    expect_false(is.null(m_perf$integrated_brier_score))

    expect_s3_class(m_perf_roc, "data.frame")

    ###################
    ### MODEL PARTS ###
    ###################

    # time dependent metrics
    m_parts_brier <- model_parts(exp, loss_function = loss_brier_score, type = "raw", output_type = "survival")
    m_parts_brier_groups <- model_parts(exp, loss_function = loss_brier_score, type = "ratio", output_type = "survival", groups = list(group1 = c("trt", "celltype"), group2 = c("age")))
    m_parts_brier_difference <- model_parts(exp, loss_function = loss_brier_score, type = "difference", output_type = "survival", variables = "age")

    plot(m_parts_brier)
    plot(m_parts_brier_groups)
    plot(m_parts_brier_difference)

    # integrated metrics
    m_parts_int_brier <- model_parts(exp, loss_function = loss_integrated_brier_score, type = "raw", output_type = "survival")
    m_parts_int_brier_groups <- model_parts(exp, loss_function = loss_integrated_brier_score, type = "ratio", output_type = "survival", groups = list(group1 = c("trt", "celltype"), group2 = c("age")))
    m_parts_int_brier_difference <- model_parts(exp, loss_function = loss_integrated_brier_score, type = "difference", output_type = "survival", variables = "age")

    plot(m_parts_int_brier)
    plot(m_parts_int_brier_groups)
    plot(m_parts_int_brier_difference)

    # c-index

    m_parts_c_ind <- model_parts(exp, loss_function = loss_one_minus_c_index, type = "raw", output_type = "risk")
    m_parts_c_ind_groups <- model_parts(exp, loss_function = loss_one_minus_c_index, type = "ratio", output_type = "risk", groups = list(group1 = c("trt", "celltype"), group2  = c("age")))
    m_parts_c_ind_difference <- model_parts(exp, loss_function = loss_one_minus_c_index, type = "difference", output_type = "risk", variables = "age")

    plot(m_parts_c_ind)
    plot(m_parts_c_ind_groups)
    plot(m_parts_c_ind_difference)

    #-- TESTS --#

    expect_s3_class(m_parts_brier, c("model_parts_survival", "surv_feature_importance", "list"))
    expect_s3_class(m_parts_brier_groups, c("model_parts_survival", "surv_feature_importance", "list"))
    expect_s3_class(m_parts_brier_difference, c("model_parts_survival", "surv_feature_importance", "list"))

    expect_s3_class(m_parts_int_brier, c("model_parts", "feature_importance_explainer", "data.frame"))
    expect_s3_class(m_parts_int_brier_groups, c("model_parts", "feature_importance_explainer", "data.frame"))
    expect_s3_class(m_parts_int_brier_difference, c("model_parts", "feature_importance_explainer", "data.frame"))

    expect_s3_class(m_parts_c_ind, c("model_parts", "feature_importance_explainer", "data.frame"))
    expect_s3_class(m_parts_c_ind_groups, c("model_parts", "feature_importance_explainer", "data.frame"))
    expect_s3_class(m_parts_c_ind_difference, c("model_parts", "feature_importance_explainer", "data.frame"))

    #####################
    ### MODEL PROFILE ###
    #####################

    m_profile <- model_profile(exp, categorical_variables = "trt")
    m_profile_subset <- model_profile(exp, variables = c("age", "diagtime"))

    # m_profile_tmp <- model_profile(exp, groups="trt")
    #
    # plot(m_profile_tmp)

    plot(m_profile)
    plot(m_profile, variable_type = "categorical")
    plot(m_profile, variables = c("diagtime", "prior"), variable_type = "numerical", numerical_plot_type = "contours", facet_ncol = 1)

    plot(m_profile_subset)

    expect_s3_class(m_profile, c("model_profile_survival", "list"))
    expect_s3_class(m_profile_subset, c("model_profile_survival", "list"))

    ###############
    ### PREDICT ###
    ###############

    sf_prediction <- predict(exp, veteran_data[1:10, ], output_type = "survival")
    chf_prediction <- predict(exp, veteran_data[1:10, ], output_type = "chf")
    risk_prediction <- predict(exp, veteran_data[1:10, ], output_type = "risk")

    #-- TESTS --#

    expect_true(is.matrix(sf_prediction))
    expect_true(is.matrix(chf_prediction))
    expect_false(is.matrix(risk_prediction))
    expect_true(is.numeric(risk_prediction))

    expect_equal(nrow(sf_prediction), 10)
    expect_equal(nrow(chf_prediction), 10)
    expect_equal(length(risk_prediction), 10)

    #####################
    ### PREDICT PARTS ###
    #####################

    p_parts_shap <- predict_parts(exp, veteran_data[1, ], type = "survshap")

    plot(p_parts_shap)

    p_parts_lime <- predict_parts(exp, veteran_data[1, ], type = "survlime")

    plot(p_parts_lime)
    plot(p_parts_lime, type = "local_importance")
    plot(p_parts_lime, show_survival_function = F)


    #-- TESTS --#

    expect_s3_class(p_parts_shap, c("predict_parts_survival", "surv_shap"))
    expect_output(print(p_parts_shap))

    expect_s3_class(p_parts_lime, c("predict_parts_survival", "surv_lime", "list"))
    expect_output(print(p_parts_lime))

    #######################
    ### PREDICT PROFILE ###
    #######################

    p_profile <- predict_profile(exp, veteran_data[1, ])
    p_profile_with_cat <- predict_profile(exp, veteran_data[1, ], categorical_variables = c("trt", "prior"))

    plot(p_profile)
    plot(p_profile_with_cat)

    #-- TESTS --#

    expect_s3_class(p_profile, c("predict_profile_survival", "surv_ceteris_paribus", "list"))
    expect_output(print(p_profile))

    expect_s3_class(p_profile_with_cat, c("predict_profile_survival", "surv_ceteris_paribus", "list"))
    expect_output(print(p_profile_with_cat))


})
