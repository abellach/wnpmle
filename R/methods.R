#' Print method for wnpmle objects
#' @param x A \code{wnpmle} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.wnpmle <- function(x, ...) {
  cat("\nWeighted NPMLE - Recurrent Events with Competing Terminal Event\n")
  cat("Type       :", x$type, "\n")
  param_name <- if (x$model == "boxcox") "rho" else "r"
  cat("Model      :", toupper(x$model), "transformation (", param_name, "=", x$rho, ")\n")
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
#' @return An object of class \code{summary.wnpmle}, printed with the
#'   coefficient table and Lambda at tau/4, tau/2 and tau.
#' @export
summary.wnpmle <- function(object, tau_grid = TRUE, ...) {
  structure(list(fit = object, tau_grid = tau_grid), class = "summary.wnpmle")
}

#' @export
print.summary.wnpmle <- function(x, ...) {
  object   <- x$fit
  tau_grid <- x$tau_grid
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
  invisible(x)
}


#' Extract coefficients from a wnpmle object
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @return A named numeric vector of regression coefficients.
#' @export
coef.wnpmle <- function(object, ...) object$coefficients


#' Extract variance-covariance matrix from a wnpmle object
#'
#' Returns the full variance-covariance matrix for (beta, Lambda).
#' To get only the beta part, use \code{vcov(fit)[1:p, 1:p]}.
#'
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @return A numeric matrix containing the variance-covariance matrix for
#'   the regression coefficients and cumulative baseline mean function.
#'   Returns an error if \code{se = "none"} was used in \code{wnpmle_fit}.
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
#' @return An object of class \code{"logLik"} with the log-likelihood value,
#'   degrees of freedom (\code{df}) equal to the number of regression
#'   coefficients, and number of observations (\code{nobs}).
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
#' @return A numeric scalar giving the AIC value.
#' @export
AIC.wnpmle <- function(object, ..., k = 2) {
  p <- length(object$coefficients)
  -2 * object$loglik + k * p
}

#' BIC for wnpmle objects
#'
#' @param object A \code{wnpmle} object.
#' @param ... Ignored.
#' @return A numeric scalar giving the BIC value.
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
#' @param tau Optional truncation time. If \code{NULL} (default), uses the
#'   maximum observed time in the data.
#' @param mark_points Logical; if \code{TRUE} (default), marks reference
#'   points at rho=1 (open circle, Ghosh-Lin) and r=1 (filled circle,
#'   proportional odds model).
#' @param file Optional path to save the plot as a PDF (e.g.
#'   \code{"loglik_profile.pdf"}). If \code{NULL} (default), plots to the
#'   current device.
#' @param width,height Size of the PDF in inches (default 6 x 4.5); only used
#'   when \code{file} is given.
#' @param verbose Logical; print progress (default \code{TRUE}).
#' @param ... Additional arguments passed to \code{\link{wnpmle_fit}}.
#'
#' @return A data frame with columns \code{model}, \code{param} and
#'   \code{loglik}, invisibly.
#'
#' @examples
#' \donttest{
#' bdata <- bladder_prep()
#' bdata_clean <- bdata[, c("id", "time", "status", "treat", "num", "size")]
#' plot_loglik(Surv(time, status) ~ treat + num + size,
#'             data = bdata_clean, id = "id")
#' }
#' @export
plot_loglik <- function(formula, data, id = "id",
                        rho_grid    = seq(0.01, 1.2, by = 0.01),
                        r_grid      = seq(0.01, 1.2, by = 0.01),
                        tau         = NULL,
                        mark_points = TRUE,
                        file        = NULL,
                        width       = 6,
                        height      = 4.5,
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
                 tau = tau,
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
                 tau = tau,
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

  i_r1_val   <- which.min(abs(r_grid   - 1))
  i_rho1_val <- which.min(abs(rho_grid - 1))

  ll_BC.new    <- c(ll_log[i_r1_val],  ll_BC)
  ll_log.new   <- c(ll_BC[i_rho1_val], ll_log)
  rho_grid.new <- c(0, rho_grid)
  r_grid.new   <- c(0, r_grid)

  # ---- plot ----
  if (!is.null(file)) pdf(file, width = width, height = height, useDingbats = FALSE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ylim_all <- range(c(ll_log.new, ll_BC.new), finite = TRUE)
  max_r    <- ceiling(max(r_grid.new)   / 0.4) * 0.4
  max_rho  <- ceiling(max(rho_grid.new) / 0.4) * 0.4
  xlim_all <- c(-max_r, max_rho)

  par(mar = c(5, 4.5, 3, 3), mgp = c(2.5, 0.5, 0), tcl = -0.25)

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
  mtext("r",             side = 1, at = -0.5 * max_r,   line = 2)
  mtext(expression(rho), side = 1, at =  0.5 * max_rho, line = 2)

  mtext("Logarithmic", side = 3, at = -0.5 * max_r,   cex = 0.9, line = 0.5)
  mtext("Box-Cox",     side = 3, at =  0.5 * max_rho, cex = 0.9, line = 0.5)

  if (mark_points) {
    i_r1   <- which.min(abs(r_grid.new   - 1))
    i_rho1 <- which.min(abs(rho_grid.new - 1))

    # labels beside the points: "r = 1" to the right, "rho = 1" above,
    # so that the curves do not run through them
    points(-r_grid.new[i_r1],    ll_log.new[i_r1],  pch = 16, cex = 1)
    text(-r_grid.new[i_r1], ll_log.new[i_r1], labels = "r = 1",
         pos = 4, offset = 0.7, xpd = NA)
    points(rho_grid.new[i_rho1], ll_BC.new[i_rho1], pch = 1,  cex = 1.5)
    text(rho_grid.new[i_rho1], ll_BC.new[i_rho1], labels = expression(rho == 1),
         pos = 3, offset = 0.9, xpd = NA)
  }

  if (!is.null(file)) {
    dev.off()
    cat("Plot saved to", file, "\n")
  }

  # ---- report optima ----
  best_rho    <- rho_grid[which.max(ll_BC)]
  best_r      <- r_grid[which.max(ll_log)]
  best_ll_bc  <- round(max(ll_BC,  na.rm = TRUE), 4)
  best_ll_log <- round(max(ll_log, na.rm = TRUE), 4)

  cat("\nOptimal transformation parameters:\n")
  tab_opt <- data.frame(
    Model   = c("Box-Cox (Ghosh-Lin at rho=1)", "Log (Prop. odds at r=1)"),
    Parameter = c(paste("rho =", best_rho), paste("r =", best_r)),
    LogLik  = c(best_ll_bc, best_ll_log)
  )
  print(tab_opt, row.names = FALSE)

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
#'   \item{lower}{Lower pointwise confidence band for Lambda.}
#'   \item{upper}{Upper pointwise confidence band for Lambda.}
#'
#' @examples
#' \donttest{
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


#' Plot method for wnpmle objects
#'
#' Plots the estimated cumulative baseline mean function Lambda(t) with
#' optional pointwise confidence bands.
#'
#' @param x A \code{wnpmle} object.
#' @param conf_bands Logical; if \code{TRUE} (default), adds pointwise 95\%
#'   confidence bands when available.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return No return value, called for side effects (produces a plot).
#' @export
plot.wnpmle <- function(x, conf_bands = TRUE, ...) {
  bl <- baseline(x)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  plot(bl$time, bl$Lambda, type = "s",
       xlab = "Time", ylab = expression(hat(Lambda)(t)),
       main = paste0("Cumulative baseline mean (", toupper(x$model),
                     ", ", if (x$model == "boxcox") "rho" else "r",
                     " = ", x$rho, ")"),
       ...)

  if (conf_bands && !is.null(bl$lower)) {
    lines(bl$time, bl$lower, lty = 2, col = "grey50")
    lines(bl$time, bl$upper, lty = 2, col = "grey50")
  }

  invisible(x)
}


#' Predict marginal mean for new covariate values
#'
#' Evaluates the estimated marginal mean number of recurrent events,
#' \eqn{\hat\mu(t \mid z) = G\{e^{\hat\beta^T z} \hat\Lambda(t)\}}, at new
#' covariate values, with pointwise confidence limits.
#'
#' @param object A \code{wnpmle} object.
#' @param newdata A data frame with the same covariates used in fitting
#'   (factors may be given as factors or character strings).
#'   If \code{NULL}, returns the estimated Lambda(t) for the baseline
#'   (all covariates = 0).
#' @param times Time points at which to evaluate the marginal mean.
#'   If \code{NULL}, uses the observed recurrent event times.
#' @param conf_int Logical; if \code{TRUE} (default) and standard errors are
#'   available, adds pointwise confidence limits.
#' @param conf_level Confidence level (default 0.95).
#' @param ... Ignored.
#'
#' @details The confidence limits are computed with the delta method from the
#'   estimated covariance matrix of \eqn{(\beta, \Lambda)} (\code{vcov(object)}),
#'   on the log scale so that they stay positive.
#'
#' @return A data frame with column \code{time} and, for each row \code{i} of
#'   \code{newdata}, columns \code{mu_i} and (if \code{conf_int}) \code{lower_i}
#'   and \code{upper_i}; or columns \code{mu} (and \code{lower}, \code{upper})
#'   for the baseline.
#' @export
predict.wnpmle <- function(object, newdata = NULL, times = NULL,
                           conf_int = TRUE, conf_level = 0.95, ...) {
  Lambda <- object$Lambda
  t_obs  <- object$event_times

  if (is.null(times)) times <- t_obs

  # interpolate Lambda at requested times (step function)
  Lambda_t <- stats::stepfun(t_obs, c(0, Lambda))(times)
  k_t      <- findInterval(times, t_obs)          # index of Lambda used (0 = before first event)
  zq       <- qnorm(1 - (1 - conf_level) / 2)
  have_se  <- conf_int && !is.null(object$vcov) && !anyNA(object$vcov)
  numcov   <- length(object$coefficients)

  log_ci <- function(mu, se) {
    ratio <- ifelse(mu > 0, se / mu, 0)
    list(lower = mu * exp(-zq * ratio), upper = mu * exp(zq * ratio))
  }

  if (is.null(newdata)) {
    out <- data.frame(time = times, mu = Lambda_t)
    if (have_se) {
      se <- ifelse(k_t > 0, sqrt(pmax(diag(object$vcov)[numcov + pmax(k_t, 1)], 0)), 0)
      ci <- log_ci(Lambda_t, se)
      out$lower <- ci$lower
      out$upper <- ci$upper
    }
    return(out)
  }

  if (!is.null(object$.terms)) {
    # same coding as in the fit (factor levels, contrasts)
    mf <- stats::model.frame(object$.terms, newdata, xlev = object$.xlevels)
    cov_mat <- model.matrix(object$.terms, mf, xlev = object$.xlevels)[, -1, drop = FALSE]
  } else {
    # objects fitted with wnpmle <= 0.1.2
    cov_mat <- model.matrix(
      stats::as.formula(paste("~", paste(object$.covars, collapse = "+"))),
      data = newdata
    )[, -1, drop = FALSE]
  }

  beta <- object$coefficients
  rho  <- object$rho
  G  <- function(x) {
    if (object$model == "boxcox") {
      if (abs(rho) < 1e-10) log(1 + x) else ((1 + x)^rho - 1) / rho
    } else log(1 + rho * x) / rho
  }
  dG <- function(x) {
    if (object$model == "boxcox") (1 + x)^(rho - 1) else 1 / (1 + rho * x)
  }

  out <- data.frame(time = times)
  for (i in seq_len(nrow(cov_mat))) {
    zi   <- cov_mat[i, ]
    e    <- exp(as.numeric(zi %*% beta))
    x    <- e * Lambda_t
    mu_t <- G(x)
    out[[paste0("mu_", i)]] <- mu_t

    if (have_se) {
      V  <- object$vcov
      gp <- dG(x)
      se <- numeric(length(times))
      for (j in which(k_t > 0)) {
        idx <- c(seq_len(numcov), numcov + k_t[j])
        g   <- c(gp[j] * x[j] * zi, gp[j] * e)   # d mu / d(beta, Lambda(t))
        se[j] <- sqrt(max(as.numeric(t(g) %*% V[idx, idx, drop = FALSE] %*% g), 0))
      }
      ci <- log_ci(mu_t, se)
      out[[paste0("lower_", i)]] <- ci$lower
      out[[paste0("upper_", i)]] <- ci$upper
    }
  }
  out
}
