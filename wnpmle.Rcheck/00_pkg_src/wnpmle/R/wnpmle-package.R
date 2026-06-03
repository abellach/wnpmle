#' wnpmle: Weighted NPMLE for Recurrent Events with a Competing Terminal Event
#'
#' Implements the weighted nonparametric maximum likelihood estimator (wNPMLE)
#' for the marginal mean function of recurrent events in the presence of a
#' competing terminal event. Two semiparametric transformation models are
#' provided: the Box-Cox transformation model and the logarithmic
#' transformation model. Both use automatic differentiation via TMB for
#' efficient and exact gradient computation.
#'
#' The main user-facing function is \code{\link{wnpmle_fit}}.
#'
#' @references
#' Bellach, A. and Kosorok, M.R. (2026). Weighted NPMLE for the marginal mean
#' of recurrent events with a competing terminal event. \emph{JASA}, to appear.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
