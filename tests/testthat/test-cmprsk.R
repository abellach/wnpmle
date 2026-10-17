# Competing risks (type = "cmprsk").
# Reference values: the stand-alone code verified on 16 Oct 2026
# (test_step3.R / test_step4.R, same data), printed to 4-5 digits.

sim_fg_data <- function() {
  # Anna's simulation code from sim.FGwnpmle.R (set.seed(666), n = 200)
  set.seed(666)
  n <- 200; p <- 0.3
  beta1 <- c(1, -0.2); beta2 <- c(0.5, 0.5)
  z1 <- rnorm(n); z2 <- rnorm(n)
  cova1 <- cbind(z1, z2)
  U <- runif(n)
  a1 <- b <- a2 <- rep(0, n); eps <- rep(0, 1)
  for (i in 1:n) {
    a1[i] <- exp(beta1 %*% cova1[i, ])
    a2[i] <- exp(beta2 %*% cova1[i, ])
    b[i]  <- 1 - (1 - p)^a1[i]
    eps[i] <- sample(x = c(1, 2), size = 1, replace = TRUE,
                     prob = c(1 - (1 - p)^a1[i], (1 - p)^a1[i]))
  }
  etime <- ifelse(eps == 1, -log(1 - (1 / p) * (1 - exp(log(1 - b * U) / a1))),
                  -(1 / a2) * log(1 - U))
  C <- runif(n, min = 0.5, max = 1)
  data.frame(id = 1:n, time = ifelse(C < etime, C, etime),
             status = ifelse(etime < C, eps, 0), z1 = z1, z2 = z2)
}

check_fit <- function(fit, est, se_f, se_s, se_a, data, family, par) {
  tg <- fit$t_grid
  expect_true(all(abs(c(fit$coefficients, fit$Lambda[tg]) - est) < 1e-4))
  expect_true(all(abs(c(fit$se, fit$se_Lambda[tg]) - se_a) < 1e-4))
  f_fisher <- wnpmle_fit(Surv(time, status) ~ z1 + z2, data = data, type = "cmprsk",
                         model = family, rho = par, tau = 0.95, se = "fisher")
  expect_true(all(abs(c(f_fisher$se, f_fisher$se_Lambda[tg]) - se_f) < 1e-4))
  f_sand <- wnpmle_fit(Surv(time, status) ~ z1 + z2, data = data, type = "cmprsk",
                       model = family, rho = par, tau = 0.95, se = "sandwich")
  expect_true(all(abs(c(f_sand$se, f_sand$se_Lambda[tg]) - se_s) < 1e-4))
}

test_that("cmprsk Box-Cox rho = 1 reproduces the verified stand-alone code", {
  skip_if_not_installed("TMB")
  d <- sim_fg_data()
  fit <- wnpmle_fit(Surv(time, status) ~ z1 + z2, data = d, type = "cmprsk",
                    model = "boxcox", rho = 1, tau = 0.95)
  expect_s3_class(fit, "wnpmle_cmprsk")
  check_fit(fit,
            est  = c(0.6773, -0.1583, 0.1125, 0.1646, 0.2563),
            se_f = c(0.15308, 0.16762, 0.02482, 0.03154, 0.05127),
            se_s = c(0.15380, 0.20588, 0.02383, 0.03067, 0.04834),
            se_a = c(0.15379, 0.20588, 0.02383, 0.03067, 0.04830),
            data = d, family = "boxcox", par = 1)
})

test_that("cmprsk log r = 1 reproduces the verified stand-alone code", {
  skip_if_not_installed("TMB")
  d <- sim_fg_data()
  fit <- wnpmle_fit(Surv(time, status) ~ z1 + z2, data = d, type = "cmprsk",
                    model = "log", rho = 1, tau = 0.95)
  check_fit(fit,
            est  = c(0.8006, -0.2155, 0.1165, 0.1791, 0.2988),
            se_f = c(0.18944, 0.18939, 0.02776, 0.03766, 0.06903),
            se_s = c(0.19393, 0.23141, 0.02688, 0.03670, 0.06446),
            se_a = c(0.19393, 0.23141, 0.02688, 0.03670, 0.06440),
            data = d, family = "log", par = 1)
})

test_that("cmprsk methods work and predict gives cumulative incidences", {
  skip_if_not_installed("TMB")
  d <- sim_fg_data()
  fit <- wnpmle_fit(Surv(time, status) ~ z1 + z2, data = d, type = "cmprsk",
                    model = "boxcox", rho = 2, tau = 0.95)
  expect_output(print(fit), "Competing Risks")
  expect_output(print(summary(fit)), "A\\(tau\\)")
  expect_equal(length(coef(fit)), 2)
  pr <- predict(fit, newdata = data.frame(z1 = c(-1, 1), z2 = 0), times = c(0.3, 0.6, 0.9))
  expect_true(all(c("cif_1", "lower_1", "upper_1", "cif_2") %in% names(pr)))
  expect_true(all(pr$cif_1 > 0 & pr$cif_1 < 1))
  expect_true(all(diff(pr$cif_1) >= 0))
  expect_true(all(pr$lower_1 <= pr$cif_1 & pr$cif_1 <= pr$upper_1))
  expect_true(all(pr$cif_2 > pr$cif_1))          # beta_1 > 0
  b0 <- predict(fit, times = c(0.3, 0.9))
  expect_true(all(b0$cif > 0 & b0$cif < 1))
  expect_error(wnpmle_fit(Surv(time, status) ~ z1, data = rbind(d, d), type = "cmprsk"),
               "one row per subject")
})
