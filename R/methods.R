#' Print method for wnpmle objects
#' @param x A \code{wnpmle} object.
#' @param ... Ignored.
#' @export
print.wnpmle <- function(x, ...) {
  cat("\nWeighted NPMLE — Recurrent Events with Competing Terminal Event\n")
  cat("Type       :", x$type, "\n")
  cat("Model      :", toupper(x$model), "transformation (rho/r =", x$rho, ")\n")
  cat("Subjects   :", x$n, "\n")
  cat("Events     : recurrent =", x$n_events["recurrent"],
      "  terminal =", x$n_events["terminal"],
      "  censored =", x$n_events["censored"], "\n")
  cat("Log-lik    :", round(x$loglik, 4), "\n")
  cat("Convergence:", x$convergence, "\n\n")
  cat("Coefficients:\n")
  tab <- cbind(Estimate = round(x$coefficients, 4),
               SE       = round(x$se, 4))
  if (!anyNA(x$se)) {
    z   <- x$coefficients / x$se
    tab <- cbind(tab,
                 "z value" = round(z, 3),
                 "Pr(>|z|)" = signif(2 * pnorm(-abs(z)), 4))
  }
  print(tab)
  invisible(x)
}


#' Summary method for wnpmle objects
#'
#' @param object A \code{wnpmle} object.
#' @param tau_grid Logical; if \code{TRUE} (default), also show Lambda at
#'   tau/4, tau/2, and tau.
#' @param ... Ignored.
#' @export
summary.wnpmle <- function(object, tau_grid = TRUE, ...) {
  print(object)

  if (tau_grid) {
    tg <- object$t_grid
    if (all(tg > 0)) {
      cat("\nCumulative baseline mean at time grid:\n")
      Lvals <- object$Lambda[tg]
      SEvals <- object$se_Lambda[tg]
      tnames <- c(
        paste0("A(tau/4) = ", round(object$tau / 4, 1)),
        paste0("A(tau/2) = ", round(object$tau / 2, 1)),
        paste0("A(tau)   = ", round(object$tau, 1))
      )
      tab <- cbind(
        Lambda = round(Lvals, 4),
        SE     = round(SEvals, 4)
      )
      rownames(tab) <- tnames
      print(tab)
    }
  }
  invisible(object)
}


#' Extract coefficients from a wnpmle object
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @export
coef.wnpmle <- function(object, ...) object$coefficients


#' Extract variance-covariance matrix from a wnpmle object
#'
#' Returns the full variance-covariance matrix for (beta, Lambda).
#' To get only the beta part, use \code{vcov(fit)[1:p, 1:p]}.
#'
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @export
vcov.wnpmle <- function(object, ...) {
  if (is.null(object$vcov))
    stop("No variance-covariance matrix available (se = 'none').")
  object$vcov
}


#' Log-likelihood for wnpmle objects
#'
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @export
logLik.wnpmle <- function(object, ...) {
  val <- object$loglik
  attr(val, "df")   <- length(object$coefficients)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' AIC for wnpmle objects
#'
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @param k Penalty per parameter (default 2 for AIC).
#' @export
AIC.wnpmle <- function(object, ..., k = 2) {
  p <- length(object$coefficients)
  -2 * object$loglik + k * p
}

#' BIC for wnpmle objects
#'
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @export
BIC.wnpmle <- function(object, ...) {
  p <- length(object$coefficients)
  -2 * object$loglik + log(object$n) * p
}


#' Log-likelihood profile plot for the transformation parameter
#'
#' Fits the model over a fine grid of transformation parameter values and
#' plots the profile log-likelihood for both the Box-Cox and logarithmic
#' transformation models on a single plot. The log transformation parameter
#' r is shown on the left (negative axis) and the Box-Cox parameter rho on
#' the right (positive axis), meeting at zero where both models coincide.
#'
#' @param formula A formula as passed to \code{\link{wnpmle_fit}}.
#' @param data A data frame.
#' @param id Name of the subject identifier column.
#' @param rho_grid A numeric vector of rho values for the Box-Cox model
#'   (default: \code{seq(0.01, 1.2, by = 0.01)}).
#' @param r_grid A numeric vector of r values for the log model
#'   (default: \code{seq(0.01, 1.2, by = 0.01)}).
#' @param mark_points Logical; if \code{TRUE} (default), marks reference
#'   points at r=0, r=1, rho=0, rho=1.
#' @param file Optional path to save the plot as a PDF (e.g.
#'   \code{"loglik_profile.pdf"}). If \code{NULL} (default), plots to the
#'   current device.
#' @param verbose Logical; print progress (default \code{TRUE}).
#' @param ... Additional arguments passed to \code{\link{wnpmle_fit}}.
#'
#' @return A data frame with columns \code{model}, \code{param} and
#'   \code{loglik}, invisibly.
#'
#' @examples
#' \dontrun{
#' bdata <- bladder_prep()
#' bdata_clean <- bdata[, c("id", "time", "status", "treat", "num", "size")]
#' plot_loglik(Surv(time, status) ~ treat + num + size,
#'             data = bdata_clean, id = "id")
#' }
#' @export
plot_loglik <- function(formula, data, id = "id",
                        rho_grid    = seq(0.01, 1.2, by = 0.01),
                        r_grid      = seq(0.01, 1.2, by = 0.01),
                        mark_points = TRUE,
                        file        = NULL,
                        verbose     = TRUE,
                        ...) {

  # ---- Box-Cox grid ----
  if (verbose) cat("Fitting Box-Cox grid (", length(rho_grid), "models)...\n")
  ll_BC     <- numeric(length(rho_grid))
  init_beta <- NULL

  for (k in seq_along(rho_grid)) {
    fit_k <- tryCatch(
      wnpmle_fit(formula, data = data, id = id,
                 model = "boxcox", rho = rho_grid[k],
                 se = "none",
                 init_beta = init_beta, ...),
      error = function(e) NULL
    )
    if (!is.null(fit_k)) {
      ll_BC[k]  <- fit_k$loglik
      init_beta <- fit_k$coefficients  # warm start
    } else {
      ll_BC[k] <- NA_real_
    }
    if (verbose && k %% 20 == 0)
      cat("  BC:", k, "/", length(rho_grid), "\n")
  }

  # ---- Log grid ----
  if (verbose) cat("Fitting log grid (", length(r_grid), "models)...\n")
  ll_log    <- numeric(length(r_grid))
  init_beta <- NULL

  for (k in seq_along(r_grid)) {
    fit_k <- tryCatch(
      wnpmle_fit(formula, data = data, id = id,
                 model = "log", rho = r_grid[k],
                 se = "none",
                 init_beta = init_beta, ...),
      error = function(e) NULL
    )
    if (!is.null(fit_k)) {
      ll_log[k]  <- fit_k$loglik
      init_beta  <- fit_k$coefficients
    } else {
      ll_log[k] <- NA_real_
    }
    if (verbose && k %% 20 == 0)
      cat("  Log:", k, "/", length(r_grid), "\n")
  }

  # match paper: use the value at param=1 of the OTHER model as starting point
  # BC curve starts at loglik of log model at r=1
  # log curve starts at loglik of BC model at rho=1
  i_r1_val   <- which.min(abs(r_grid   - 1))
  i_rho1_val <- which.min(abs(rho_grid - 1))

  ll_BC.new    <- c(ll_log[i_r1_val],  ll_BC)
  ll_log.new   <- c(ll_BC[i_rho1_val], ll_log)
  rho_grid.new <- c(0, rho_grid)
  r_grid.new   <- c(0, r_grid)

  # ---- plot ----
  if (!is.null(file)) pdf(file, width = 3.5, height = 3.5, useDingbats = FALSE)

  ylim_all <- range(c(ll_log.new, ll_BC.new), finite = TRUE)
  max_r    <- ceiling(max(r_grid.new)   / 0.4) * 0.4
  max_rho  <- ceiling(max(rho_grid.new) / 0.4) * 0.4
  xlim_all <- c(-max_r, max_rho)

  par(mar = c(6, 5, 4, 2), mgp = c(1.25, 0.22, 0), tcl = -0.18)

  plot(NA, xlim = xlim_all, ylim = ylim_all,
       xlab = "", ylab = "Log-likelihood", axes = FALSE)
  box()

  lines(-r_grid.new,   ll_log.new, lwd = 2, lty = 2)
  lines(rho_grid.new,  ll_BC.new,  lwd = 2)
  abline(v = 0, lty = 3, col = "grey60")

  ticks_left  <- round(seq(-max_r,  0,       by = 0.4), 2)
  ticks_right <- round(seq(0,       max_rho, by = 0.4), 2)
  axis(1, at     = c(ticks_left, ticks_right[-1]),
          labels = c(abs(ticks_left), ticks_right[-1]))
  axis(2)

  mtext("Transformation parameter", side = 1, line = 3.5)
  mtext("r",             side = 1, at = -0.7 * max_r,   line = 1.8)
  mtext(expression(rho), side = 1, at =  0.7 * max_rho, line = 1.8)

  if (mark_points) {
    i_r1   <- which.min(abs(r_grid.new   - 1))
    i_rho1 <- which.min(abs(rho_grid.new - 1))

    points(-r_grid.new[i_r1],    ll_log.new[i_r1],  pch = 16, cex = 1)
    text(-r_grid.new[i_r1] + 0.14, ll_log.new[i_r1], labels = "r = 1", adj = 0)

    points(rho_grid.new[i_rho1], ll_BC.new[i_rho1], pch = 1,  cex = 1)
    text(rho_grid.new[i_rho1] + 0.14, ll_BC.new[i_rho1],
         labels = expression(rho == 1), adj = 0)
  }

  if (!is.null(file)) {
    dev.off()
    cat("Plot saved to", file, "\n")
  }

  # ---- report optima ----
  best_rho <- rho_grid[which.max(ll_BC)]
  best_r   <- r_grid[which.max(ll_log)]
  cat("Optimal rho (Box-Cox):", best_rho,
      "  loglik:", round(max(ll_BC, na.rm = TRUE), 4), "\n")
  cat("Optimal r   (Log)    :", best_r,
      "  loglik:", round(max(ll_log, na.rm = TRUE), 4), "\n")

  invisible(data.frame(
    model  = c(rep("log", length(r_grid)),   rep("boxcox", length(rho_grid))),
    param  = c(r_grid,                        rho_grid),
    loglik = c(ll_log,                        ll_BC)
  ))
}


#' Extract the estimated baseline mean function
#'
#' Returns a data frame with the estimated cumulative baseline mean function
#' Lambda(t) and its increments lambda(t), together with standard errors and
#' pointwise 95% confidence intervals.
#'
#' @param object A \code{wnpmle} object.
#' @param conf_level Confidence level for the pointwise intervals (default 0.95).
#' @param ... Ignored.
#'
#' @return A data frame with columns:
#'   \item{time}{Recurrent event times.}
#'   \item{lambda}{Estimated baseline increments.}
#'   \item{Lambda}{Estimated cumulative baseline mean.}
#'   \item{se_Lambda}{Standard error of Lambda (if SE was estimated).}
#'   \item{lower}{Lower confidence band for Lambda.}
#'   \item{upper}{Upper confidence band for Lambda.}
#'
#' @examples
#' \dontrun{
#' bdata <- bladder_prep()
#' bdata_clean <- bdata[, c("id", "time", "status", "treat", "num", "size")]
#' fit <- wnpmle_fit(Surv(time, status) ~ treat + num + size,
#'                   data = bdata_clean, id = "id", model = "log", rho = 1)
#' bl <- baseline(fit)
#' head(bl)
#' }
#' @export
baseline <- function(object, conf_level = 0.95, ...) {
  if (!inherits(object, "wnpmle"))
    stop("object must be of class 'wnpmle'")

  z    <- qnorm(1 - (1 - conf_level) / 2)
  out  <- data.frame(
    time   = object$event_times,
    lambda = object$lambda,
    Lambda = object$Lambda
  )

  if (!anyNA(object$se_Lambda)) {
    out$se_Lambda <- object$se_Lambda
    out$lower     <- pmax(object$Lambda - z * object$se_Lambda, 0)
    out$upper     <- object$Lambda + z * object$se_Lambda
  }

  out
}


#' Predict marginal mean for new covariate values
#'
#' @param object A \code{wnpmle} object.
#' @param newdata A data frame with the same covariates used in fitting.
#'   If \code{NULL}, returns the estimated Lambda(t) for the baseline
#'   (all covariates = 0).
#' @param times Time points at which to evaluate the marginal mean.
#'   If \code{NULL}, uses the observed recurrent event times.
#' @param ... Ignored.
#' @return A data frame with columns \code{time} and one column per row of
#'   \code{newdata} (or a single column \code{mu} for baseline).
#' @export
predict.wnpmle <- function(object, newdata = NULL, times = NULL, ...) {
  Lambda <- object$Lambda
  t_obs  <- object$event_times

  if (is.null(times)) times <- t_obs

  # interpolate Lambda at requested times (step function)
  Lambda_t <- stats::stepfun(t_obs, c(0, Lambda))(times)

  if (is.null(newdata)) {
    return(data.frame(time = times, mu = Lambda_t))
  }

  cov_mat <- model.matrix(
    stats::as.formula(paste("~", paste(object$.covars, collapse = "+"))),
    data = newdata
  )[, -1, drop = FALSE]

  beta <- object$coefficients
  out  <- data.frame(time = times)

  for (i in seq_len(nrow(cov_mat))) {
    eta <- as.numeric(cov_mat[i, ] %*% beta)
    if (object$model == "boxcox") {
      rho <- object$rho
      if (abs(rho) < 1e-10) {
        mu_t <- exp(eta) * Lambda_t
      } else {
        mu_t <- ((1 + exp(eta) * Lambda_t)^rho - 1) / rho
      }
    } else {
      r    <- object$rho
      mu_t <- log(1 + r * exp(eta) * Lambda_t) / r
    }
    out[[paste0("mu_", i)]] <- mu_t
  }
  out
}
