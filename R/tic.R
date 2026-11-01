#' Takeuchi information criterion for the weighted NPMLE
#'
#' Computes Takeuchi's information criterion (TIC) for a fitted \code{wnpmle}
#' model. The weighted NPMLE maximizes an inverse probability of censoring
#' weighted log-likelihood, for which the information equality need not hold.
#' TIC therefore replaces the AIC penalty \eqn{2p} by
#' \eqn{2\,\mathrm{tr}(\mathcal{I}_n^{-1}\widehat\Sigma_n)}, where
#' \eqn{\mathcal{I}_n} is the observed information and \eqn{\widehat\Sigma_n}
#' the censoring-corrected meat of the sandwich estimator:
#' \deqn{\mathrm{TIC} = -2\,\ell(\hat\theta) + 2\,\mathrm{tr}(\mathcal{I}_n^{-1}\widehat\Sigma_n).}
#'
#' The trace is computed from the full covariance matrices of
#' \eqn{(\beta, \Lambda)}:
#' \eqn{\mathrm{tr}(\mathcal{I}_n^{-1}\widehat\Sigma_n) = \mathrm{tr}(V_S V_F^{-1})},
#' with \eqn{V_S = \mathcal{I}_n^{-1}\widehat\Sigma_n\mathcal{I}_n^{-1}} the
#' sandwich and \eqn{V_F = \mathcal{I}_n^{-1}} the inverse Fisher information.
#' The trace is invariant to the parametrization of the baseline (jumps or
#' cumulative values).
#'
#' @param object A \code{wnpmle} object fitted with a sandwich variance
#'   estimator (\code{se = "sandwich_adj"}, the default, or \code{"sandwich"}).
#' @param fisher Optional \code{wnpmle} object for the same model and data,
#'   fitted with \code{se = "fisher"}. If \code{NULL} (default), the model is
#'   refitted with \code{se = "fisher"}; this requires that the data used in
#'   the original call are still available.
#' @param details Logical; if \code{TRUE}, also return the penalty, the number
#'   of parameters, AIC and the log-likelihood (default \code{FALSE}).
#' @param ... Ignored.
#'
#' @return A numeric scalar giving the TIC value. With \code{details = TRUE},
#'   a named vector with \code{TIC}, \code{penalty} (the trace
#'   \eqn{\mathrm{tr}(\mathcal{I}_n^{-1}\widehat\Sigma_n)}), \code{p} (the
#'   number of parameters, including the baseline jumps), \code{AIC}
#'   (\eqn{-2\ell + 2p} with this \code{p}) and \code{loglik}.
#'
#' @details
#' For a correctly specified model with a genuine likelihood the penalty
#' equals \code{p} asymptotically, and TIC reduces to AIC. Since all Box-Cox
#' and logarithmic transformation models have the same number of parameters,
#' AIC compares log-likelihoods only, whereas TIC can give different
#' transformations different penalties. The trace is estimated, so TIC is
#' noisier than AIC in small samples.
#'
#' The penalty should not exceed \code{p} by much; a penalty clearly larger
#' than \code{p} points to a problem with the variance estimate.
#'
#' @references
#' Takeuchi, K. (1976). Distribution of informational statistics and a
#' criterion of model fitting. \emph{Suri-Kagaku (Mathematical Sciences)},
#' 153, 12-18 (in Japanese).
#'
#' Konishi, S. and Kitagawa, G. (1996). Generalised information criteria in
#' model selection. \emph{Biometrika}, 83, 875-890.
#'
#' @seealso \code{\link{wnpmle_fit}}, \code{\link{AIC.wnpmle}}
#'
#' @examples
#' \dontrun{
#'   dat <- bladder_prep()
#'   fit <- wnpmle_fit(Surv(time, status) ~ treat + num + size, data = dat,
#'                     id = "id", model = "boxcox", rho = 1)
#'   TIC(fit)
#'   TIC(fit, details = TRUE)
#' }
#'
#' @export
TIC <- function(object, ...) UseMethod("TIC")

#' @rdname TIC
#' @export
TIC.wnpmle <- function(object, fisher = NULL, details = FALSE, ...) {
  se_used <- object$call$se
  if (is.null(se_used)) se_used <- "sandwich_adj"
  if (!se_used %in% c("sandwich_adj", "sandwich"))
    stop("'object' must be fitted with se = \"sandwich_adj\" or \"sandwich\".")

  if (is.null(fisher)) {
    cl <- object$call
    cl$se <- "fisher"
    fisher <- eval(cl, parent.frame())
  }

  Vs <- object$vcov
  Vf <- fisher$vcov
  if (is.null(Vs) || is.null(Vf))
    stop("Full covariance matrices are required in both fits.")
  if (!identical(dim(Vs), dim(Vf)))
    stop("Covariance matrices of the two fits have different dimensions.")
  if (abs(object$loglik - fisher$loglik) > 1e-6 * max(1, abs(object$loglik)))
    warning("The two fits have different log-likelihoods; are model and data the same?")

  penalty <- tryCatch(sum(diag(Vs %*% solve(Vf))),
                      error = function(e) {
                        warning("Inverse Fisher covariance not invertible: ",
                                conditionMessage(e))
                        NA_real_
                      })
  p   <- nrow(Vs)
  tic <- -2 * object$loglik + 2 * penalty
  if (is.finite(penalty) && penalty > 1.2 * p)
    warning(sprintf("TIC penalty (%.1f) clearly exceeds p (%d); check the variance estimate.",
                    penalty, p))
  if (details)
    return(c(TIC = tic, penalty = penalty, p = p,
             AIC = -2 * object$loglik + 2 * p, loglik = object$loglik))
  tic
}
