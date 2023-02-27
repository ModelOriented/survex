test_that("CD/AUC works", {

    rotterdam <- survival::rotterdam
    rotterdam$year <- NULL
    rotterdam <- rotterdam[sample(1:nrow(rotterdam), 300), ]

    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )

    coxph_explainer <- explain(cox_rotterdam_rec, verbose = FALSE)

    times <- coxph_explainer$times
    surv <- coxph_explainer$predict_survival_function(coxph_explainer$model, coxph_explainer$data, times)
    auc <- loss_one_minus_cd_auc(y_true = coxph_explainer$y, surv = surv, times = times)

    expect_equal(length(auc), length(times))
    expect_lt(sum(is.na(auc)), length(times))
    expect_true(all(auc >= 0, na.rm = TRUE))
    expect_true(all(auc <= 1, na.rm = TRUE))


})


test_that("C-index works", {
    rotterdam <- survival::rotterdam
    rotterdam$year <- NULL

    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )

    coxph_explainer <- explain(cox_rotterdam_rec, y = survival::Surv(rotterdam$rtime, rotterdam$recur), verbose = FALSE)

    times <- coxph_explainer$times
    risk <- coxph_explainer$predict_function(coxph_explainer$model, coxph_explainer$data)
    c_ind <- loss_one_minus_c_index(y_true = coxph_explainer$y, risk = risk, times = times)

    expect_lte(c_ind, 1)
    expect_gte(c_ind, 0)

})


test_that("Brier score works", {
    rotterdam <- survival::rotterdam
    rotterdam$year <- NULL

    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )

    coxph_explainer <- explain(cox_rotterdam_rec, y = survival::Surv(rotterdam$rtime, rotterdam$recur), verbose = FALSE)


    times <- coxph_explainer$times
    surv <- coxph_explainer$predict_survival_function(coxph_explainer$model, coxph_explainer$data, times)
    bs <- brier_score(y_true = coxph_explainer$y, surv = surv, times = times)

    expect_equal(length(bs), length(times))
    expect_lt(sum(is.na(bs)), length(times))
    expect_true(all(bs >= 0, na.rm = TRUE))
    expect_true(all(bs <= 1, na.rm = TRUE))

})

test_that("integration of metrics works", {
    rotterdam <- survival::rotterdam
    rotterdam$year <- NULL

    cox_rotterdam_rec <-
        survival::coxph(
            survival::Surv(rtime, recur) ~ .,
            data = rotterdam[, !colnames(rotterdam) %in% c("year", "dtime", "death")],
            model = TRUE,
            x = TRUE,
            y = TRUE
        )

    coxph_explainer <- explain(cox_rotterdam_rec, y = survival::Surv(rotterdam$rtime, rotterdam$recur), verbose = FALSE)


    times <- coxph_explainer$times
    surv <- coxph_explainer$predict_survival_function(coxph_explainer$model, coxph_explainer$data, times)

    ibs <- loss_integrate(loss_brier_score, normalization = "t_max")
    bs <- ibs(y_true = coxph_explainer$y, surv = surv, times = times)

    expect_length(bs, 1)
    expect_lte(bs, 1)
    expect_gte(bs, 0)

    # ibs <- loss_integrate(loss_brier_score, normalization = "survival")
    # bs <- ibs(y_true = coxph_explainer$y, surv = surv, times = times)
    #
    # expect_length(bs, 1)
    # expect_lte(bs, 1)
    # expect_gte(bs, 0)

})
