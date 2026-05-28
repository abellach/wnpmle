test_that("bladder_prep returns correct structure", {
  skip_if_not_installed("survival")
  bdata <- bladder_prep()

  expect_s3_class(bdata, "data.frame")
  expect_true(all(c("id", "time", "status", "rx", "size", "number") %in% names(bdata)))
  expect_true(all(bdata$status %in% 0:2))
  expect_true(sum(bdata$status == 1) > 0)
  expect_true(sum(bdata$status == 2) > 0)
})

test_that("wnpmle_fit returns wnpmle object with log model", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  bdata <- bladder_prep()
  fit <- wnpmle_fit(
    survival::Surv(time, status) ~ rx + size + number,
    data  = bdata,
    id    = "id",
    model = "log",
    rho   = 1,
    se    = "sandwich"
  )

  expect_s3_class(fit, "wnpmle")
  expect_named(fit$coefficients, c("rx", "size", "number"))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$convergence, "relative convergence (4)")
  expect_equal(length(fit$Lambda), fit$n_events["recurrent"])
})

test_that("wnpmle_fit returns wnpmle object with boxcox model", {
  skip_if_not_installed("survival")
  skip_if_not_installed("TMB")

  bdata <- bladder_prep()
  fit <- wnpmle_fit(
    survival::Surv(time, status) ~ rx + size + number,
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
    survival::Surv(time, status) ~ rx + size + number,
    data  = bdata,
    id    = "id",
    model = "log",
    rho   = 1,
    se    = "none"
  )

  pred0 <- predict(fit)
  expect_s3_class(pred0, "data.frame")
  expect_true(all(pred0$mu >= 0))

  newdat <- data.frame(rx = c(1, 2), size = c(1, 2), number = c(1, 3))
  pred1 <- predict(fit, newdata = newdat)
  expect_equal(ncol(pred1), 3)  # time + mu_1 + mu_2
})
