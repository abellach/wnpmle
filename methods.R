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
  attr(val, "df") <- length(object$coefficients)
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
#' Fits the model over a grid of transformation parameter values (rho or r)
#' and plots the profile log-likelihood. Useful for selecting the optimal
#' transformation.
#'
#' @param formula A formula as passed to \code{\link{wnpmle_fit}}.
#' @param data A data frame.
#' @param id Name of the subject identifier column.
#' @param model Transformation model: \code{"boxcox"} or \code{"log"}.
#' @param rho_grid A numeric vector of rho/r values to evaluate
#'   (default: \code{seq(0.5, 5, by = 0.5)} for Box-Cox,
#'   \code{seq(0.5, 3, by = 0.5)} for log).
#' @param se Variance estimation for each fit (default \code{"none"} for speed).
#' @param plot Logical; if \code{TRUE} (default), plot the profile.
#' @param ... Additional arguments passed to \code{\link{wnpmle_fit}}.
#'
#' @return A data frame with columns \code{rho} and \code{loglik}, invisibly.
#'
#' @examples
#' \dontrun{
#' bdata <- bladder_prep()
#' bdata_clean <- bdata[, c("id", "time", "status", "treat", "num", "size")]
#' plot_loglik(Surv(time, status) ~ treat + num + size,
#'             data = bdata_clean, id = "id", model = "boxcox")
#' }
#' @export
plot_loglik <- function(formula, data, id = "id",
                        model     = c("boxcox", "log"),
                        rho_grid  = NULL,
                        se        = "none",
                        plot      = TRUE,
                        ...) {
  model <- match.arg(model)

  if (is.null(rho_grid)) {
    rho_grid <- if (model == "boxcox") seq(0.5, 5, by = 0.5) else
                                        seq(0.5, 3, by = 0.5)
  }

  logliks <- numeric(length(rho_grid))
  cat("Fitting", length(rho_grid), "models...\n")

  for (i in seq_along(rho_grid)) {
    cat("  rho =", rho_grid[i], "\n")
    fit_i <- tryCatch(
      wnpmle_fit(formula, data = data, id = id,
                 model = model, rho = rho_grid[i], se = se, ...),
      error = function(e) NULL
    )
    logliks[i] <- if (!is.null(fit_i)) fit_i$loglik else NA_real_
  }

  result <- data.frame(rho = rho_grid, loglik = logliks)

  if (plot) {
    xlab <- if (model == "boxcox") expression(rho) else "r"
    plot(rho_grid, logliks, type = "b", pch = 19,
         xlab = xlab, ylab = "Log-likelihood",
         main = paste("Profile log-likelihood —",
                      ifelse(model == "boxcox", "Box-Cox", "Log"), "model"))
    best <- rho_grid[which.max(logliks)]
    abline(v = best, lty = 2, col = "red")
    legend("topright",
           legend = paste0("Best ", ifelse(model == "boxcox", "rho", "r"),
                           " = ", best),
           lty = 2, col = "red", bty = "n")
  }

  invisible(result)
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
