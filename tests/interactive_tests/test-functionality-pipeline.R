
test_that("Manual testing pipeline", {
skip("For manual testing only")
veteran <- survival::veteran

cph <- survival::coxph(survival::Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE, y = TRUE)

rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ .,
                                    data = veteran,
                                    respect.unordered.factors = TRUE,
                                    num.trees = 100,
                                    mtry = 3,
                                    max.depth = 5)

rsf_src <- randomForestSRC::rfsrc(Surv(time, status) ~ .,
                                        data = veteran)

cph_exp <- explain(cph, verbose = FALSE)
rsf_ranger_exp <- explain(rsf_ranger, data = veteran[, -c(3, 4)], y = Surv(veteran$time, veteran$status), verbose = FALSE)
rsf_src_exp <- explain(rsf_src, verbose = FALSE) # data=veteran[,-c(3,4)], y=Surv(veteran$time, veteran$status))

#########################
### model_performance ###
#########################

# metrics vizualization
cph_model_performance <- model_performance(cph_exp)
plot(cph_model_performance)

rsf_ranger_model_performance <- model_performance(rsf_ranger_exp)
rsf_src_model_performance <- model_performance(rsf_src_exp)

plot(rsf_ranger_model_performance, cph_model_performance, rsf_src_model_performance, metrics_type = "scalar")
plot(rsf_ranger_model_performance, cph_model_performance, rsf_src_model_performance, metrics_type = "scalar", colors = c("red", "green", "blue"))
plot(rsf_ranger_model_performance, cph_model_performance, rsf_src_model_performance)

# selecting only a subset of metrics
plot(rsf_src_model_performance, rsf_ranger_model_performance, cph_model_performance,
     metrics_type = "time_dependent", metrics = c("Brier score"))

# ROC curves
cph_model_performance_roc <- model_performance(cph_exp, type = "roc", times = c(100, 500, 1200))
ranger_model_performance_roc <- model_performance(rsf_ranger_exp, type = "roc", times = c(100, 500, 1200))
src_model_performance_roc <- model_performance(rsf_src_exp, type = "roc", times = c(100, 500, 1200))

plot(cph_model_performance_roc, ranger_model_performance_roc, src_model_performance_roc)


###################
### model_parts ###
###################

# time-dependent metrics
cph_model_parts_brier <- model_parts(cph_exp)
plot(cph_model_parts_brier)

rsf_ranger_model_parts <- model_parts(rsf_ranger_exp)
plot(cph_model_parts_brier, rsf_ranger_model_parts)

# specifying loss function brier
rsf_src_model_parts_brier <- model_parts(rsf_src_exp, loss_function = loss_brier_score,
                                   output_type = "survival")
plot(rsf_src_model_parts_brier)

# c/d auc
rsf_src_model_parts_auc <- model_parts(rsf_src_exp, loss_function = loss_one_minus_cd_auc,
                                      output_type = "survival", B = 1)
plot(rsf_src_model_parts_auc, max_vars = 3)


# integrated AUC
cph_model_parts_int_auc <- model_parts(cph_exp, loss_function = loss_one_minus_integrated_cd_auc, B = 2)
plot(cph_model_parts_int_auc)

# integrated brier
cph_model_parts_int_brier <- model_parts(cph_exp, loss_function = loss_integrated_brier_score, B = 10)
plot(cph_model_parts_int_brier)

ranger_model_parts_int_brier <- model_parts(rsf_ranger_exp, loss_function = loss_integrated_brier_score, B = 10)
src_model_parts_int_brier <- model_parts(rsf_src_exp, loss_function = loss_integrated_brier_score, B = 10)
plot(cph_model_parts_int_brier, ranger_model_parts_int_brier, src_model_parts_int_brier)

# C index - default DALEX
cph_model_parts_dalex <- model_parts(cph_exp, loss_function = loss_one_minus_c_index,
                               output_type = "risk")
class(cph_model_parts_dalex)
plot(cph_model_parts_dalex)

#####################
### model_profile ###
#####################

# numerical variables
cph_model_profile <- model_profile(cph_exp, output_type = "survival")
plot(cph_model_profile)
plot(cph_model_profile, colors = c("red", "green", "blue"))

rsf_ranger_model_profile <- model_profile(rsf_ranger_exp, output_type = "survival")
plot(rsf_ranger_model_profile)

# categorical variables
cph_model_profile_cat <- model_profile(cph_exp, output_type = "survival", variables = "celltype")
plot(cph_model_profile_cat, variable_type = "categorical")

###############
### predict ###
###############

predict(cph_exp, newdata = veteran[1, ], output_type = "risk")
predict(cph_exp, newdata = veteran[2, ], output_type = "survival")
predict(cph_exp, newdata = veteran[3, ], output_type = "chf")

predict(rsf_ranger_exp, newdata = veteran[1, ], output_type = "risk")
predict(rsf_ranger_exp, newdata = veteran[2, ], output_type = "survival")
predict(rsf_ranger_exp, newdata = veteran[3, ], output_type = "chf")

#####################
### predict_parts ###
#####################

# SHAP
cph_predict_parts_survshap <- predict_parts(cph_exp, new_observation = veteran[1, -c(3, 4)])
plot(cph_predict_parts_survshap)
plot(cph_predict_parts_survshap, colors = c("cyan", "pink", "brown", "gray", "green", "yellow"))

rsf_ranger_predict_parts_survshap <- predict_parts(rsf_ranger_exp, new_observation = veteran[1, -c(3, 4)])
plot(cph_predict_parts_survshap, rsf_ranger_predict_parts_survshap)


# LIME
cph_predict_parts_survlime <- predict_parts(cph_exp, new_observation = veteran[1, -c(3, 4)], type = "survlime")
plot(cph_predict_parts_survlime)
plot(cph_predict_parts_survlime, type = "local_importance")
plot(cph_predict_parts_survlime, show_survival_function = FALSE)

#######################
### predict_profile ###
#######################

cph_predict_profile <- predict_profile(cph_exp, veteran[2, -c(3, 4)])
plot(cph_predict_profile, colors = c("#ff0000", "#00ff00", "#0000ff"))
plot(cph_predict_profile)


cph_predict_profile_cat <- predict_profile(cph_exp, veteran[2, -c(3, 4)], variables = c("celltype"))
plot(cph_predict_profile_cat, variable_type = "categorical", colors = c("#ff0000", "#00ff00", "#0000ff"))
plot(cph_predict_profile_cat, variable_type = "categorical")

})
