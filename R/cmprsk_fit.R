# ============================================================================
# Competing risks: wnpmle_fit(..., type = "cmprsk")
#
# Weighted NPMLE for the subdistribution of a competing risk
# (Bellach, Kosorok, Rueschendorf and Fine, 2019), Box-Cox and logarithmic
# transformation models. Built from the stand-alone code verified on
# 16 Oct 2026 (cmprsk_fit.R, test_step12/3/4.R):
#   estimates and inverse Fisher   src/TMB/fn_cmprsk_BC.hpp,    fn_cmprsk_log.hpp
#   subject-wise scores            src/TMB/score_cmprsk_BC.hpp, score_cmprsk_log.hpp
#   censoring correction psi       src/TMB/corr_cmprsk_BC.hpp,  corr_cmprsk_log.hpp
# The recurrent-events code (wnpmle_fit.R, helpers.R) is not used here.
# ============================================================================

# TMB model name for competing risks, as listed in src/wnpmle.cpp
.wnpmle_tmb_model_cmprsk <- function(model, part = c("fn", "score", "corr")) {
  part <- match.arg(part)
  paste0(part, "_cmprsk_", if (model == "boxcox") "BC" else "log")
}

# Data preparation (as cmprsk_prep() in the verified stand-alone code).
# One row per subject; status 0 = censored, 1 = event of interest,
# 2 = competing event. Rows are ordered by time, events before censorings at
# ties. Subjects observed beyond tau, and the last subjects once the
# Kaplan-Meier estimate of the censoring distribution drops below 0.05, are
# censored just after tau (as in Bellach et al., 2019).
.cmprsk_prep <- function(time, status, Z, tau) {
  Z <- as.matrix(Z); n <- length(time)
  o <- order(time, status == 0)
  time <- time[o]; status <- status[o]; Z <- Z[o, , drop = FALSE]

  # Kaplan-Meier of the censoring distribution, by position
  kmc <- numeric(n)
  kmc[1] <- (1 - 1 / n)^(status[1] == 0)
  if (n > 1) for (i in 2:n) kmc[i] <- kmc[i - 1] * (1 - 1 / (n - i + 1))^(status[i] == 0)

  nbig <- sum(kmc < .05)
  if (nbig > 0) {
    idx <- (n - nbig + 1):n
    time[idx] <- tau + 1e-5 * seq_len(nbig); status[idx] <- 0; kmc[idx] <- kmc[n - nbig]
  }
  nbig2 <- sum(time > tau)
  if (nbig2 > 0) {
    idx <- (n - nbig2 + 1):n
    time[idx] <- tau + 1e-5 * seq_len(nbig2); status[idx] <- 0; kmc[idx] <- kmc[n - nbig2]
  }

  ind  <- seq_len(n)
  ind1 <- ind[status == 1]; ind2 <- ind[status == 2]; indc <- ind[status == 0]
  num1 <- length(ind1); num2 <- length(ind2)
  if (num1 == 0) stop("No events of interest (status == 1) at or before tau.")

  w <- matrix(0, num2, num1)
  for (j in seq_len(num2)) w[j, ] <- (ind1 > ind2[j]) * (kmc[ind1] / kmc[ind2[j]])

  # starting values for the jumps (as in the simulation scripts)
  compr <- numeric(num1)
  compr[1] <- sum(status[1:ind1[1]])
  if (num1 > 1) for (j in 2:num1) compr[j] <- sum(status[ind1[j - 1]:ind1[j]] == 2)
  compr  <- compr / kmc[ind1]
  salpha <- log(1 / ((n - ind1) + cumsum(compr) * kmc[ind1]))

  j2 <- rep(-1L, n); j2[ind2] <- seq_len(num2) - 1L
  list(n = n, time = time, status = status, Z = Z, kmc = kmc, order = o,
       ind1 = ind1, ind2 = ind2, indc = indc, num1 = num1, num2 = num2,
       cov1 = Z[ind1, , drop = FALSE], cov2 = Z[ind2, , drop = FALSE], w = w,
       idx_all = findInterval(ind, ind1), salpha = salpha, j2 = j2)
}

.wnpmle_fit_cmprsk <- function(formula, data, id, model, rho, tau, se,
                               init_beta, control, silent, cl) {
  # ---- 1. Parse formula and data ----
  resp <- as.character(formula[[2]])
  if (!grepl("Surv", resp[1]) || length(resp) < 3)
    stop("Response must be a Surv object: Surv(time, status)")
  time_var <- resp[2]; status_var <- resp[3]
  if (!time_var %in% names(data))
    stop("Time variable '", time_var, "' not found in data.")
  if (!status_var %in% names(data))
    stop("Status variable '", status_var, "' not found in data.")
  time_col   <- as.numeric(data[[time_var]])
  status_col <- as.integer(data[[status_var]])
  if (!all(status_col %in% 0:2))
    stop("status must be 0 (censored), 1 (event of interest), or 2 (competing event).")
  if (!is.null(id) && id %in% names(data) && anyDuplicated(data[[id]]))
    stop("type = \"cmprsk\" needs one row per subject, but column '", id,
         "' has repeated values.")

  rhs_formula <- formula
  rhs_formula[[2]] <- NULL
  cov_mat   <- model.matrix(rhs_formula, data)[, -1, drop = FALSE]
  rhs_terms <- stats::delete.response(stats::terms(rhs_formula, data = data))
  rhs_xlev  <- stats::.getXlevels(rhs_terms, stats::model.frame(rhs_terms, data))
  covars    <- colnames(cov_mat)
  numcov    <- ncol(cov_mat)

  if (is.null(tau)) tau <- max(time_col)
  prep <- .cmprsk_prep(time_col, status_col, cov_mat, tau)
  num1 <- prep$num1

  # ---- 2. Fit (estimates; inverse Fisher in (beta, lambda) space) ----
  dat <- list(model = .wnpmle_tmb_model_cmprsk(model, "fn"),
              cova = prep$Z, idx_all = as.integer(prep$idx_all),
              cov1 = prep$cov1, cov2 = prep$cov2, w = prep$w)
  if (model == "boxcox") dat$rho <- as.numeric(rho) else dat$r <- as.numeric(rho)
  beta_init  <- if (!is.null(init_beta)) init_beta else rep(0, numcov)
  parameters <- list(beta = beta_init, alpha = prep$salpha)

  obj <- TMB::MakeADFun(dat, parameters, DLL = "wnpmle", silent = silent)
  opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control = control)
  if (opt$convergence != 0) warning("nlminb: ", opt$message)

  beta_hat   <- opt$par[seq_len(numcov)]
  lambda_hat <- exp(opt$par[numcov + seq_len(num1)])
  Lambda_hat <- cumsum(lambda_hat)
  names(beta_hat) <- covars
  loglik <- -opt$objective

  time1 <- prep$time[prep$ind1]
  t_grid <- c(tau4 = sum(time1 <= tau / 4), tau2 = sum(time1 <= tau / 2),
              tau = sum(time1 <= tau))

  # ---- 3. Variance: (beta, lambda) -> (beta, Lambda) ----
  vcov_mat  <- NULL
  se_beta   <- rep(NA_real_, numcov)
  se_Lambda <- rep(NA_real_, num1)
  if (se != "none") {
    H <- obj$he(opt$par)
    V <- tryCatch(chol2inv(chol(H)), error = function(e) {
      warning("Hessian is not positive definite; using pseudoinverse.")
      MASS::ginv(H)
    })
    J      <- diag(c(rep(1, numcov), lambda_hat))
    breadi <- J %*% V %*% J
    Mtrafo <- rbind(cbind(diag(numcov), matrix(0, numcov, num1)),
                    cbind(matrix(0, num1, numcov), outer(seq_len(num1), seq_len(num1), ">=") * 1))

    if (se == "fisher") {
      vcov_mat <- Mtrafo %*% breadi %*% t(Mtrafo)
    } else {
      # subject-wise scores of the negative log-likelihood w.r.t. (beta, lambda)
      sdat <- list(model = .wnpmle_tmb_model_cmprsk(model, "score"),
                   cova = prep$Z, idx_all = as.integer(prep$idx_all),
                   status = as.integer(prep$status), j2 = as.integer(prep$j2), w = prep$w)
      if (model == "boxcox") sdat$rho <- as.numeric(rho) else sdat$r <- as.numeric(rho)
      sobj  <- TMB::MakeADFun(sdat, list(beta = beta_hat, alpha = log(lambda_hat)),
                              DLL = "wnpmle", ADreport = TRUE, silent = silent)
      gradi <- sobj$gr(c(beta_hat, log(lambda_hat)))
      gradi[, numcov + seq_len(num1)] <-
        sweep(gradi[, numcov + seq_len(num1), drop = FALSE], 2, lambda_hat, "/")

      if (se == "sandwich_adj") {
        # censoring correction; gradi is on the negative log-likelihood scale,
        # so the corrected meat uses gradi - psi (= -(eta_i + kappa_i))
        cdat <- list(model = .wnpmle_tmb_model_cmprsk(model, "corr"),
                     cova = prep$Z, ind1 = as.integer(prep$ind1 - 1L),
                     ind2 = as.integer(prep$ind2 - 1L), indc = as.integer(prep$indc - 1L),
                     w = prep$w)
        if (model == "boxcox") cdat$rho <- as.numeric(rho) else cdat$r <- as.numeric(rho)
        cobj  <- TMB::MakeADFun(cdat, list(beta = beta_hat, alpha = log(lambda_hat)),
                                DLL = "wnpmle", silent = silent)
        gradi <- gradi - cobj$report(c(beta_hat, log(lambda_hat)))$psi
      }
      vcov_mat <- Mtrafo %*% breadi %*% crossprod(gradi) %*% breadi %*% t(Mtrafo)
    }
    se_all    <- sqrt(pmax(diag(vcov_mat), 0))
    se_beta   <- se_all[seq_len(numcov)]
    se_Lambda <- se_all[numcov + seq_len(num1)]
    names(se_beta) <- covars
  }

  structure(
    list(
      coefficients = beta_hat,
      se           = se_beta,
      Lambda       = Lambda_hat,
      se_Lambda    = se_Lambda,
      lambda       = lambda_hat,
      event_times  = time1,
      loglik       = loglik,
      vcov         = vcov_mat,
      model        = model,
      type         = "cmprsk",
      se_type      = se,
      rho          = rho,
      tau          = tau,
      n            = prep$n,
      n_events     = c(event = num1, competing = prep$num2,
                       censored = length(prep$indc)),
      zeng_lin     = FALSE,
      t_grid       = t_grid,
      convergence  = opt$message,
      call         = cl,
      .covars      = covars,
      .terms       = rhs_terms,
      .xlevels     = rhs_xlev
    ),
    class = c("wnpmle_cmprsk", "wnpmle")
  )
}


# ---- Methods for competing-risks fits ----------------------------------------

#' @rdname print.wnpmle
#' @export
print.wnpmle_cmprsk <- function(x, ...) {
  cat("\nWeighted NPMLE - Competing Risks (subdistribution of event type 1)\n")
  param_name <- if (x$model == "boxcox") "rho" else "r"
  cat("Model      :", toupper(x$model), "transformation (", param_name, "=", x$rho, ")\n")
  cat("Subjects   :", x$n, "\n")
  cat("Events     : event of interest =", x$n_events["event"],
      "  competing =", x$n_events["competing"],
      "  censored =", x$n_events["censored"], "\n")
  cat("Log-lik    :", round(x$loglik, 4), "\n")
  cat("Convergence:", x$convergence, "\n\n")
  cat("Coefficients:\n")
  tab <- cbind(Estimate = round(x$coefficients, 4), SE = round(x$se, 4))
  if (!anyNA(x$se)) {
    z   <- x$coefficients / x$se
    tab <- cbind(tab, "z value" = round(z, 3), "Pr(>|z|)" = signif(2 * pnorm(-abs(z)), 4))
  }
  print(tab)
  invisible(x)
}

#' @rdname summary.wnpmle
#' @export
summary.wnpmle_cmprsk <- function(object, tau_grid = TRUE, ...) {
  structure(list(fit = object, tau_grid = tau_grid), class = "summary.wnpmle_cmprsk")
}

#' @export
print.summary.wnpmle_cmprsk <- function(x, ...) {
  object <- x$fit
  print(object)
  if (x$tau_grid) {
    tg <- object$t_grid
    if (all(tg > 0)) {
      cat("\nCumulative baseline function at time grid:\n")
      tab <- cbind(Lambda = round(object$Lambda[tg], 4), SE = round(object$se_Lambda[tg], 4))
      rownames(tab) <- c(paste0("A(tau/4) = ", round(object$tau / 4, 2)),
                         paste0("A(tau/2) = ", round(object$tau / 2, 2)),
                         paste0("A(tau)   = ", round(object$tau, 2)))
      print(tab)
    }
  }
  invisible(x)
}

#' @rdname plot.wnpmle
#' @export
plot.wnpmle_cmprsk <- function(x, conf_bands = TRUE, ...) {
  bl <- baseline(x)
  plot(bl$time, bl$Lambda, type = "s",
       xlab = "Time", ylab = expression(hat(Lambda)(t)),
       main = paste0("Cumulative baseline function (", toupper(x$model), ", ",
                     if (x$model == "boxcox") "rho" else "r", " = ", x$rho, ")"),
       ...)
  if (conf_bands && !is.null(bl$lower)) {
    lines(bl$time, bl$lower, lty = 2, col = "grey50")
    lines(bl$time, bl$upper, lty = 2, col = "grey50")
  }
  invisible(x)
}

#' @rdname predict.wnpmle
#' @details For competing-risks fits (\code{type = "cmprsk"}), \code{predict}
#'   returns the cumulative incidence of the event of interest,
#'   \eqn{\hat F_1(t \mid z) = 1 - \exp[-G\{e^{\hat\beta^T z} \hat\Lambda(t)\}]},
#'   in columns \code{cif_i}; the confidence limits are those of
#'   \eqn{G\{e^{\beta^T z}\Lambda(t)\}} (delta method, log scale) mapped through
#'   \eqn{1 - \exp(-\cdot)}, so they stay within (0, 1).
#' @export
predict.wnpmle_cmprsk <- function(object, newdata = NULL, times = NULL,
                                  conf_int = TRUE, conf_level = 0.95, ...) {
  # the marginal-mean prediction gives G(e^{beta'z} Lambda(t)) and its limits;
  # the cumulative incidence is 1 - exp(-G(.)), a monotone map
  p <- NextMethod()
  rho <- object$rho
  G <- function(x) {
    if (object$model == "boxcox") {
      if (abs(rho) < 1e-10) log(1 + x) else ((1 + x)^rho - 1) / rho
    } else log(1 + rho * x) / rho
  }
  # without newdata the marginal-mean method returns Lambda(t) itself (z = 0),
  # so G has to be applied first; with newdata it already returns G(.)
  to_cif <- if (is.null(newdata)) function(v) 1 - exp(-G(v)) else function(v) 1 - exp(-v)
  nm <- names(p)
  for (k in setdiff(nm, "time")) p[[k]] <- to_cif(p[[k]])
  names(p) <- sub("^mu", "cif", nm)
  p
}
