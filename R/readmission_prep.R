#' Prepare the colorectal cancer readmission data for wnpmle analysis
#'
#' Reads the hospital readmission data of 403 colorectal cancer patients
#' (Gonzalez et al., 2005), shipped with \pkg{wnpmle} as a copy of the
#' \code{readmission} data from the \pkg{frailtypack} package, and converts
#' them to the format used by \code{\link{wnpmle_fit}}: one row per
#' readmission (status 1) and one final row per patient, which is either death
#' (status 2) or censoring (status 0).
#'
#' @param tau Optional end of follow-up (days). If given, readmissions after
#'   \code{tau} are removed and patients still under observation at
#'   \code{tau} (including deaths after \code{tau}) are censored at \code{tau}.
#'   Default \code{NULL}: full follow-up.
#'
#' @return A data frame with columns:
#'   \item{id}{Patient identifier.}
#'   \item{time}{Days since surgery of the readmission, death or censoring.}
#'   \item{status}{0 = censored, 1 = readmission (recurrent event),
#'     2 = death (terminal event).}
#'   \item{chemo}{Chemotherapy: factor with levels \code{"NonTreated"},
#'     \code{"Treated"}.}
#'   \item{sex}{Factor with levels \code{"Male"}, \code{"Female"}.}
#'   \item{dukes}{Dukes' tumour stage: factor with levels \code{"A-B"},
#'     \code{"C"}, \code{"D"}.}
#'   \item{charlson}{Charlson comorbidity index (time-dependent): factor with
#'     levels \code{"0"}, \code{"1-2"}, \code{"3"}.}
#'
#' @details
#' The data contain 861 rows for 403 patients: 458 readmissions, 109 deaths
#' and 294 censored follow-ups. The Charlson index may change between
#' readmissions; each row carries the value recorded for that interval.
#'
#' @source
#' Copy of the \code{readmission} data from the \pkg{frailtypack} package
#' (Rondeau et al.; GPL (>= 2)).
#'
#' @references
#' Gonzalez, J.R., Fernandez, E., Moreno, V., Ribes, J., Peris, M., Navarro, M.,
#' Cambray, M. and Borras, J.M. (2005). Sex differences in hospital readmission
#' among colorectal cancer patients. \emph{Journal of Epidemiology and
#' Community Health}, 59(6), 506-511. \doi{10.1136/jech.2004.028902}
#'
#' @examples
#' rdata <- readmission_prep()
#' head(rdata)
#' table(rdata$status)
#' table(readmission_prep(tau = 1460)$status)
#'
#' @export
readmission_prep <- function(tau = NULL) {
  file <- system.file("extdata", "readmission.csv", package = "wnpmle",
                      mustWork = TRUE)
  raw <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)

  covs <- c("chemo", "sex", "dukes", "charlson")

  # one row per readmission
  rec <- raw[raw$event == 1, , drop = FALSE]
  rec <- data.frame(id = rec$id, time = rec[["t.stop"]], status = 1L,
                    rec[, covs], stringsAsFactors = FALSE)

  # last row of each patient: death (2) or censoring (0)
  last <- raw[!duplicated(raw$id, fromLast = TRUE), , drop = FALSE]
  last <- data.frame(id = last$id, time = last[["t.stop"]],
                     status = ifelse(last$death == 1, 2L, 0L),
                     last[, covs], stringsAsFactors = FALSE)

  if (!is.null(tau)) {
    # end of follow-up at tau: drop later readmissions, censor at tau
    rec  <- rec[rec$time <= tau, , drop = FALSE]
    late <- last$time > tau
    last$time[late]   <- tau
    last$status[late] <- 0L
  }

  out <- rbind(rec, last)
  out <- out[order(out$id, out$time, -out$status), ]
  rownames(out) <- NULL

  out$chemo    <- factor(out$chemo,    levels = c("NonTreated", "Treated"))
  out$sex      <- factor(out$sex,      levels = c("Male", "Female"))
  out$dukes    <- factor(out$dukes,    levels = c("A-B", "C", "D"))
  out$charlson <- factor(out$charlson, levels = c("0", "1-2", "3"))
  out
}
