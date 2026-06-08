#' Fit Weighted NPMLE for Survival Data with Recurrent or Competing Events
#'
#' Estimates the weighted nonparametric maximum likelihood estimator (wNPMLE)
#' for the marginal mean function of recurrent events in the presence of a
#' competing terminal event. Two transformation models are supported:
#' the Box-Cox model and the logarithmic model.
#'
#' @param formula A formula of the form \code{Surv(time, status) ~ covariates},
#'   where \code{status} takes value \code{1} for a recurrent event,
#'   \code{2} for the terminal event (e.g. death), and \code{0} for censoring.
#' @param data A data frame containing the variables in \code{formula} and
#'   an \code{id} column identifying subjects.
#' @param id Name of the subject identifier column in \code{data}
#'   (default: \code{"id"}).
#' @param type Type of analysis. Currently only \code{"recurrent"} is
#'   supported (recurrent events with a competing terminal event).
#'   Future versions will add \code{"cmprsk"} (competing risks)
#'   and \code{"cmprsk_ltrc"} (competing risks with left truncation
#'   and right censoring).
#' @param model Transformation model: \code{"boxcox"} (default) or
#'   \code{"log"}.
#' @param rho Transformation parameter. For \code{model = "boxcox"}, this is
#'   the Box-Cox parameter rho (default 1, i.e. linear/proportional means
#'   model). For \code{model = "log"}, this is the parameter r in
#'   \code{G(x) = log(1 + r*x) / r} (default 1).
#' @param tau Follow-up truncation time. Kaplan-Meier censoring weights are
#'   truncated at \code{tau}. If \code{NULL} (default), uses the maximum
#'   observed time.
#' @param se Variance estimation method: \code{"sandwich_adj"} (default,
#'   sandwich with censoring correction), \code{"sandwich_adjold"} (same
#'   correction using original loop-based code, for validation only),
#'   \code{"sandwich"} (plain sandwich), \code{"fisher"}, or \code{"none"}.
#' @param init_beta Initial values for regression coefficients (default: all
#'   zeros).
#' @param control A list of control parameters passed to \code{\link[stats]{nlminb}}.
#' @param silent Suppress TMB output (default: \code{TRUE}).
#'
#' @return An object of class \code{"wnpmle"} with components:
#'   \item{coefficients}{Estimated regression coefficients (beta).}
#'   \item{se}{Standard errors for the regression coefficients.}
#'   \item{Lambda}{Estimated cumulative baseline mean function, evaluated at
#'     the recurrent event times.}
#'   \item{lambda}{Estimated baseline increments.}
#'   \item{event_times}{Recurrent event times at which Lambda is estimated.}
#'   \item{loglik}{Log-likelihood at the optimum.}
#'   \item{vcov}{Full variance-covariance matrix for (beta, Lambda).}
#'   \item{model}{Transformation model used.}
#'   \item{rho}{Transformation parameter used.}
#'   \item{tau}{Truncation time used.}
#'   \item{n}{Number of subjects.}
#'   \item{n_events}{Named vector: recurrent events, terminal events, censored.}
#'   \item{convergence}{Convergence message from \code{nlminb}.}
#'   \item{call}{The matched call.}
#'
#' @details
#' **Analysis types:** Currently only \code{type = "recurrent"} is implemented,
#' for the recurrent events with competing terminal event setting. Support for
#' \code{"competing_risks"} (Bellach, Kosorok, Rüschendorf and Fine, 2019) and
#' \code{"competing_risks_ltrc"} (left truncation and right censoring) will be
#' added in future versions.
#'
#' **Transformation models:** The two models differ in the link function \eqn{G}:
#' \itemize{
#'   \item **Box-Cox**: \eqn{G(x) = ((1 + x)^\rho - 1) / \rho}, with
#'     \eqn{G(x) \to \log(1+x)} as \eqn{\rho \to 0}.
#'   \item **Logarithmic**: \eqn{G(x) = \log(1 + r \cdot x) / r}, with
#'     \eqn{G(x) \to x} as \eqn{r \to 0}.
#' }
#' Both are estimated via automatic differentiation using TMB. The sandwich
#' variance estimator accounts for the estimation of the censoring weights.
#' The censoring-corrected sandwich adds an influence function correction
#' for the Kaplan-Meier censoring weight estimation.
#'
#' @references
#' Bellach, A. and Kosorok, M.R. (2026). Weighted NPMLE for the marginal mean
#' of recurrent events with a competing terminal event.
#' \url{https://arxiv.org/abs/2605.25934}.
#'
#' Bellach, A., Kosorok, M.R., Rüschendorf, L. and Fine, J.P. (2019). Weighted NPMLE for the
#' subdistribution of a competing risk. \emph{Journal of the American
#' Statistical Association}, 114(525), 259-270.
#'
#' @examples
#'  \dontrun{
#'   library(survival)
#'   data("bladder2", package = "survival")
#'   bladder2_prepped <- bladder_prep()
#' 
#'   fit <- wnpmle_fit(Surv(time, status) ~ treat + num + size,
#'                     data = bladder2_prepped, id = "id",
#'                     model = "log", rho = 1)
#'   summary(fit)
#' }
#'
#' @export
wnpmle_fit <- function(formula,
                       data,
                       id      = "id",
                       type    = c("recurrent", "competing_risks",
                                   "competing_risks_ltrc"),
                       model   = c("boxcox", "log"),
                       rho     = 1,
                       tau     = NULL,
                       se      = c("sandwich", "sandwich_corrected",
                                   "fisher", "none"),
                       init_beta = NULL,
                       control   = list(),
                       silent    = TRUE) {

  cl    <- match.call()
  type  <- match.arg(type)
  model <- match.arg(model)
  se    <- match.arg(se)

  if (type != "recurrent") {
    stop(
      "type = \"", type, "\" is not yet implemented.\n",
      "Currently only type = \"recurrent\" is supported.\n",
      "Support for competing risks will be added in a future version."
    )
  }

  # ---- 1. Parse formula and data ----
  # Extract time and status variable names from the Surv() call in the formula
  resp <- as.character(formula[[2]])
  # resp is like c("Surv", "time", "status") or c("survival::Surv", "time", "status")
  if (!grepl("Surv", resp[1]) || length(resp) < 3)
    stop("Response must be a Surv object: Surv(time, status)")

  time_var   <- resp[2]
  status_var <- resp[3]

  if (!time_var %in% names(data))
    stop("Time variable '", time_var, "' not found in data.")
  if (!status_var %in% names(data))
    stop("Status variable '", status_var, "' not found in data.")

  time_col   <- data[[time_var]]
  status_col <- as.integer(data[[status_var]])

  # build covariate matrix from RHS of formula
  rhs_formula <- formula
  rhs_formula[[2]] <- NULL
  cov_mat <- model.matrix(rhs_formula, data)[, -1, drop = FALSE]

  if (!id %in% names(data))
    stop("Column '", id, "' not found in data.")
  id_vec <- as.integer(factor(data[[id]]))

  if (!all(status_col %in% 0:2))
    stop("status must be 0 (censored), 1 (recurrent event), or 2 (terminal event).")

  # ---- 2. Build working data frame ----
  mydata <- data.frame(
    id     = id_vec,
    time   = time_col,
    status = status_col
  )
  mydata <- cbind(mydata, as.data.frame(cov_mat))
  covars <- colnames(cov_mat)

  mydata <- mydata[order(mydata$time), ]
  row.names(mydata) <- NULL
  mydata$ind <- seq_len(nrow(mydata))

  n    <- nrow(mydata)
  numi <- length(unique(mydata$id))

  mydata$status1 <- as.integer(mydata$status == 1)
  mydata$status2 <- as.integer(mydata$status == 2)
  mydata$status0 <- as.integer(mydata$status == 0)
  mydata$dimind  <- as.integer(mydata$status != 1)

  num1 <- sum(mydata$status1)
  num2 <- sum(mydata$status2)
  numc <- sum(mydata$status0)
  n02  <- num2 + numc

  if (num1 == 0) stop("No recurrent events (status == 1) found.")
  zeng_lin <- (num2 == 0)   # no terminal events: reduce to Zeng-Lin (plain NPMLE)
  if (zeng_lin && se != "none") {
    message("No terminal events found: fitting Zeng-Lin model (plain NPMLE). ",
            "Switching variance estimator to 'sandwich'.")
    se <- "sandwich"
  }

  # ---- 3. Truncation time ----
  if (is.null(tau)) tau <- max(mydata$time)

  # ---- 4. KM censoring weights ----
  # When num2==0 (Zeng-Lin): kmc=1 for all (weights are known, not estimated)
  if (zeng_lin) {
    mydata$kmc <- 1
  } else {
    risk      <- numi - cumsum(mydata$dimind[seq_len(n - 1)])
    kmc_start <- (1 - 1 / (numi - 1))^(mydata$status[1] == 0)
    mydata$kmc <- kmc_start *
      c(1, cumprod((1 - 1 / risk)^mydata$status0[seq(2, n)]))
    i_tau <- max(which(mydata$time <= tau), 0L)
    if (i_tau > 0L) mydata$kmc[i_tau:n] <- mydata$kmc[i_tau]
  }

  # ---- 5. Subsets ----
  M   <- mydata
  M1  <- M[M$status == 1, , drop = FALSE]
  M02 <- M[M$status != 1, , drop = FALSE];  M02$idM02 <- seq_len(nrow(M02))
  Mc  <- M02[M02$status == 0, , drop = FALSE]
  M2  <- M02[M02$status == 2, , drop = FALSE]

  cova   <- data.matrix(M[, covars, drop = FALSE])
  cov1   <- cova[M$status1 == 1, , drop = FALSE]
  cov2   <- cova[M$status2 == 1, , drop = FALSE]
  cov02  <- cova[M$status1 != 1, , drop = FALSE]
  covc   <- cova[M$status0 == 1, , drop = FALSE]
  numcov <- ncol(cova)

  # indicator matrices (outer product form, consistent with simulation code)
  M3 <- outer(M1$ind,  M1$ind,  FUN = ">=") * 1L
  M5 <- outer(Mc$ind,  M1$ind,  FUN = ">=") * 1L
  M6 <- outer(M2$ind,  M1$ind,  FUN = ">=") * 1L

  wnew <- matrix(0, nrow = num2, ncol = num1)
  for (i in seq_len(num2))
    wnew[i, ] <- (M1$ind > M2$ind[i]) * (M1$kmc / M2$kmc[i])

  # ---- 6. TMB index vectors ----
  ind1  <- M1$ind
  ind02 <- M02$ind
  ind2  <- M2$ind

  if (model == "boxcox") {
    # BC cpp: idx02 values -1..num1-1
    idx02 <- as.integer(findInterval(ind02, ind1)) - 1L
    idx2  <- as.integer(findInterval(ind2,  ind1)) - 1L
    stopifnot(min(idx02) >= -1L, max(idx02) <= num1 - 1L)
  } else {
    # log cpp: idx02 values 0..num1 (no -1L)
    idx02 <- as.integer(findInterval(ind02, ind1))
    idx2  <- as.integer(findInterval(ind2,  ind1)) - 1L
    stopifnot(min(idx02) >= 0L, max(idx02) <= num1)
  }
  if (length(idx2) > 0) stopifnot(min(idx2) >= -1L, max(idx2) <= num1 - 1L)

  kmc1_tmb <- as.numeric(M1$kmc)
  kmc2_tmb <- as.numeric(M2$kmc)

  # ---- 7. TMB data and parameters ----
  rho_val  <- as.numeric(rho)
  dll_name <- if (model == "boxcox") "fn_BC_tmb" else "fn_log_tmb"

  data_tmb <- list(
    cov1  = cov1,
    cov2  = cov2,
    cov02 = cov02,
    idx02 = idx02,
    idx2  = idx2,
    kmc1  = kmc1_tmb,
    kmc2  = kmc2_tmb
  )
  if (model == "boxcox") {
    data_tmb$rho <- rho_val
  } else {
    data_tmb$r <- rho_val
  }

  beta_init <- if (!is.null(init_beta)) init_beta else rep(0, numcov)
  parameters <- list(
    beta  = beta_init,
    alpha = rep(log(1 / num1), num1)
  )

  # ---- 8. Fit ----
  obj <- TMB::MakeADFun(data_tmb, parameters, DLL = dll_name, silent = silent)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

  est        <- obj$env$parList(opt$par)
  beta_hat   <- est$beta
  lambda_hat <- exp(est$alpha)
  Lambda_hat <- cumsum(lambda_hat)
  names(beta_hat) <- covars

  loglik <- -opt$objective

  # time points for Lambda
  t4 <- sum(M1$time < tau / 4)
  t2 <- sum(M1$time < tau / 2)
  t1 <- sum(M1$time < tau + 0.1)

  # transformation for Lambda (from alpha to lambda parameterisation)
  J       <- diag(c(rep(1, numcov), lambda_hat))
  Mtrafo  <- rbind(
    cbind(diag(numcov),            matrix(0, numcov, num1)),
    cbind(matrix(0, num1, numcov), M3)
  )

  # ---- 9. Variance estimation ----
  vcov_mat <- NULL
  se_beta  <- rep(NA_real_, numcov)
  se_Lambda <- rep(NA_real_, num1)

  if (se != "none") {
    H_theta <- obj$he(opt$par)
    V_theta <- tryCatch(
      chol2inv(chol(H_theta)),
      error = function(e) {
        warning("Hessian is not positive definite; using pseudoinverse.")
        MASS::ginv(H_theta)
      }
    )
    breadi   <- J %*% V_theta %*% t(J)
    breadinv <- Mtrafo %*% breadi %*% t(Mtrafo)

    if (se == "fisher") {
      vcov_mat <- breadinv
    }

    if (se %in% c("sandwich", "sandwich_adj", "sandwich_adjold")) {
      Lambda <- as.numeric(M3 %*% lambda_hat)
      Lamc   <- as.numeric(M5 %*% lambda_hat)
      Lam2   <- as.numeric(M6 %*% lambda_hat)
      beta   <- beta_hat

      if (se == "sandwich_adjold") {
        gradi <- .compute_score_old(
          model, rho_val, numcov, num1, numi, n02, num2,
          cov1, cov2, cov02, covc, beta, lambda_hat, Lambda,
          Lamc, Lam2, wnew, M1, M2, Mc
        )
      } else {
        gradi <- .compute_score(
          model, rho_val, numcov, num1, numi, n02, num2,
          cov1, cov2, cov02, covc, beta, lambda_hat, Lambda,
          Lamc, Lam2, wnew, M1, M2, Mc
        )
      }

      if (se %in% c("sandwich_adj", "sandwich_adjold") && !zeng_lin) {
        if (se == "sandwich_adjold") {
          psi_subj <- .censoring_correction_old(
            model, rho_val, numcov, num1, numi, n02, num2,
            cov2, beta, lambda_hat, Lambda, wnew,
            M1, M2, M02
          )
        } else {
          psi_subj <- .censoring_correction(
            model, rho_val, numcov, num1, numi, n02, num2,
            cov2, beta, lambda_hat, Lambda, wnew,
            M1, M2, M02
          )
        }
        gradi_eff <- gradi + psi_subj
      } else {
        gradi_eff <- gradi
      }

      grad <- array(0, dim = c(numcov + num1, numcov + num1, numi))
      for (i in seq_len(numi))
        grad[, , i] <- tcrossprod(gradi_eff[i, ])
      meat     <- apply(grad, 1:2, sum)
      sandw    <- Mtrafo %*% breadi %*% meat %*% t(breadi) %*% t(Mtrafo)
      vcov_mat <- sandw
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
      event_times  = M1$time,
      loglik       = loglik,
      vcov         = vcov_mat,
      model        = model,
      type         = type,
      rho          = rho,
      tau          = tau,
      n            = numi,
      n_events     = c(recurrent = num1, terminal = num2, censored = numc),
      zeng_lin     = zeng_lin,
      t_grid       = c(tau4 = t4, tau2 = t2, tau = t1),
      convergence  = opt$message,
      call         = cl,
      # internal objects needed for predict/plot
      .M1          = M1,
      .M02         = M02,
      .covars      = covars
    ),
    class = "wnpmle"
  )
}
