test_that("bladder_prep returns correct structure", {
  skip_if_not_installed("survival")
  bdata <- bladder_prep()

  expect_s3_class(bdata, "data.frame")
  expect_true(all(c("id", "time", "status") %in% names(bdata)))
  expect_true(all(bdata$status %in% 0:2))
  expect_true(sum(bdata$status == 1) > 0)
  expect_true(sum(bdata$status == 2) > 0)
})

test_that("wnpmle_fit returns wnpmle object with log model", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  bdata <- bladder_prep()
  fit <- wnpmle_fit(
    Surv(time, status) ~ treat + num + size,
    data  = bdata,
    id    = "id",
    model = "log",
    rho   = 1,
    se    = "sandwich"
  )

  expect_s3_class(fit, "wnpmle")
  expect_named(fit$coefficients, c("treat", "num", "size"))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$convergence, "relative convergence (4)")
  expect_equal(length(fit$Lambda), unname(fit$n_events["recurrent"]))
})

test_that("wnpmle_fit returns wnpmle object with boxcox model", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  bdata <- bladder_prep()
  fit <- wnpmle_fit(
    Surv(time, status) ~ treat + num + size,
    data  = bdata,
    id    = "id",
    model = "boxcox",
    rho   = 1,
    se    = "fisher"
  )

  expect_s3_class(fit, "wnpmle")
  expect_true(is.finite(fit$loglik))
})

test_that("predict.wnpmle works for baseline and newdata", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  bdata <- bladder_prep()
  fit <- wnpmle_fit(
    Surv(time, status) ~ treat + num + size,
    data  = bdata,
    id    = "id",
    model = "log",
    rho   = 1,
    se    = "none"
  )

  pred0 <- predict(fit)
  expect_s3_class(pred0, "data.frame")
  expect_true(all(pred0$mu >= 0))

  newdat <- data.frame(treat = c(1, 2), num = c(1, 3), size = c(1, 2))
  pred1 <- predict(fit, newdata = newdat)
  expect_equal(ncol(pred1), 3)  # time + mu_1 + mu_2
})

test_that("sandwich standard errors are finite for both models", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  bdata <- bladder_prep()
  for (m in c("boxcox", "log")) {
    for (s in c("sandwich", "sandwich_adj")) {
      fit <- wnpmle_fit(Surv(time, status) ~ treat + num + size,
                        data = bdata, id = "id", model = m, rho = 1, se = s)
      expect_true(all(is.finite(fit$se)))
      expect_true(all(fit$se > 0))
    }
  }
})

test_that("bladder results are unchanged (Box-Cox rho = 1, tau = 59)", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  # Reference values from wnpmle 0.1.3 (vignette comparison with Ghosh-Lin).
  bdata <- bladder_prep()
  fit <- wnpmle_fit(Surv(time, status) ~ treat + num + size,
                    data = bdata, id = "id", model = "boxcox", rho = 1, tau = 59)
  expect_true(all(abs(unname(fit$coefficients) - c(-0.536, 0.173, -0.005)) < 1e-3))
  expect_true(all(abs(unname(fit$se)           - c( 0.268, 0.059,  0.068)) < 1e-3))
})
