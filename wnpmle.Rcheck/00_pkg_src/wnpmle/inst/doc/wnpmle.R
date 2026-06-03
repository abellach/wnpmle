## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 5,
  fig.height = 4.5
)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("wnpmle")

## ----eval = FALSE-------------------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("abellach/wnpmle")

## ----load---------------------------------------------------------------------
library(wnpmle)

bdata       <- bladder_prep()
bdata_clean <- bdata[, c("id", "time", "status", "treat", "num", "size")]
head(bdata_clean)

## ----fit-bc-------------------------------------------------------------------
fit_bc <- wnpmle_fit(
  Surv(time, status) ~ treat + num + size,
  data  = bdata_clean,
  id    = "id",
  model = "boxcox",
  rho   = 1,
  tau   = 59,
  se    = "sandwich"
)
summary(fit_bc)
bl_bc <- baseline(fit_bc)
plot(bl_bc$time, bl_bc$Lambda, type = "s", lwd = 2,
     xlab = "Time (months)", ylab = expression(hat(Lambda)(t)),
     main = "Cumulative baseline mean (Ghosh-Lin)")
lines(bl_bc$time, bl_bc$lower, lty = 2, col = "grey50")
lines(bl_bc$time, bl_bc$upper, lty = 2, col = "grey50")
AIC(fit_bc)
BIC(fit_bc)

## ----fit-log------------------------------------------------------------------
fit_log <- wnpmle_fit(
  Surv(time, status) ~ treat + num + size,
  data  = bdata_clean,
  id    = "id",
  model = "log",
  rho   = 1,
  tau   = 59,
  se    = "sandwich"
)
summary(fit_log)
bl_log <- baseline(fit_log)
plot(bl_log$time, bl_log$Lambda, type = "s", lwd = 2,
     xlab = "Time (months)", ylab = expression(hat(Lambda)(t)),
     main = "Cumulative baseline mean (proportional odds)")
lines(bl_log$time, bl_log$lower, lty = 2, col = "grey50")
lines(bl_log$time, bl_log$upper, lty = 2, col = "grey50")
AIC(fit_log)
BIC(fit_log)

## ----predict------------------------------------------------------------------
newdat <- data.frame(treat = c(0, 1), num = c(1, 1), size = c(1, 1))

pred <- predict(fit_bc, newdata = newdat,
                times = seq(1, 50, by = 1))

plot(pred$time, pred$mu_1, type = "s", lwd = 2,
     xlab = "Time (months)", ylab = "Marginal mean number of recurrences",
     ylim = range(pred[, -1]))
lines(pred$time, pred$mu_2, lwd = 2, lty = 2, col = "firebrick")
legend("topleft", legend = c("Placebo", "Thiotepa"),
       lty = c(1, 2), col = c("black", "firebrick"), bty = "n")

## ----loglik-plot, eval = FALSE------------------------------------------------
# plot_loglik(
#   Surv(time, status) ~ treat + num + size,
#   data     = bdata_clean,
#   id       = "id",
#   tau      = 59,
#   rho_grid = seq(0.01, 1.2, by = 0.01),
#   r_grid   = seq(0.01, 1.2, by = 0.01)
# )

