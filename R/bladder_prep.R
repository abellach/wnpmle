#' Prepare bladder cancer data for wnpmle analysis
#'
#' Prepares the \code{bladder1} dataset from the \pkg{survival} package for
#' use with \code{\link{wnpmle_fit}}. Subjects randomised to pyridoxine are
#' excluded (only placebo and thiotepa arms are retained). Status codes are
#' recoded to 0 (censored), 1 (recurrent tumour), 2 (terminal/dropout), and
#' event times beyond \code{tau} are truncated to \code{tau}.
#'
#' @param tau Follow-up truncation time (default: 59 months). Event times
#'   beyond \code{tau} are set to \code{tau}.
#'
#' @return A data frame with columns:
#'   \item{id}{Subject identifier (integer).}
#'   \item{status}{Event status: 0 = censored, 1 = recurrent event,
#'     2 = terminal event.}
#'   \item{status0}{Indicator: 1 if censored.}
#'   \item{status1}{Indicator: 1 if recurrent event.}
#'   \item{status2}{Indicator: 1 if terminal event.}
#'   \item{time}{Event time (months).}
#'   \item{treat}{Treatment indicator: 1 = thiotepa, 0 = placebo.}
#'   \item{num}{Initial number of tumours.}
#'   \item{size}{Initial tumour size (cm).}
#'
#' @details
#' The bladder cancer trial (Veterans Administration Cooperative Urological
#' Research Group) randomised patients to placebo, pyridoxine, or thiotepa.
#' Only the placebo and thiotepa arms are used here (pyridoxine arm excluded).
#' Subjects with event times equal to \code{tau} and status 0 have their
#' times truncated to \code{tau}, consistent with administrative censoring.
#'
#' @references
#' Byar, D.P. (1980). The Veterans Administration study of chemoprophylaxis
#' for recurrent stage I bladder tumors. In \emph{Bladder Tumors and Other
#' Topics in Urological Oncology}, 363-370. Plenum, New York.
#'
#' @examples
#' bdata <- bladder_prep(tau = 59)
#' head(bdata)
#' table(bdata$status)
#'
#' @export
bladder_prep <- function(tau = 59) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for bladder_prep().")

  # keep placebo (ids 1-128) and thiotepa (ids 214-294) arms only
  data <- rbind(survival::bladder1[1:128, ],
                survival::bladder1[214:294, ])

  # re-index ids after removing pyridoxine arm
  data$id[129:209] <- data$id[129:209] - 32

  # treatment indicator: 1 = thiotepa, 0 = placebo
  data$treatment <- as.numeric(data$treatment == "thiotepa")

  # recode status 3 -> 2 (death/removal both treated as terminal)
  data$status[data$status == 3] <- 2

  numi <- data$id[nrow(data)]

  # rename columns
  my <- data
  colnames(my) <- c("id", "treat", "num", "size", "recur",
                    "start", "time", "status", "rtumor", "rsize", "enum")

  status  <- my$status
  status0 <- as.integer(status == 0)
  status1 <- as.integer(status == 1)
  status2 <- as.integer(status == 2)

  mydata <- data.frame(
    id      = my$id,
    status  = status,
    status0 = status0,
    status1 = status1,
    status2 = status2,
    time    = my$time,
    treat   = my$treat,
    num     = my$num,
    size    = my$size
  )

  # add subjects censored at tau who are missing from bladder1
  mydata2 <- data.frame(
    id      = c(13, 15, 16, 19, 24, 34, 44, 51, 72),
    status  = rep(0, 9),
    status0 = rep(1, 9),
    status1 = rep(0, 9),
    status2 = rep(0, 9),
    time    = rep(64, 9),
    treat   = c(rep(0, 7), 1, 1),
    num     = c(1, 2, 1, 1, 1, 5, 3, 8, 3),
    size    = c(1, 3, 1, 4, 6, 1, 1, 1, 1)
  )
  mydata <- rbind(mydata, mydata2)

  # terminal/censored indicator
  mydata$dimind <- as.integer(mydata$status != 1)

  # order by id and time within id
  mydata <- mydata[order(mydata$id, mydata$time), ]
  rownames(mydata) <- seq_len(nrow(mydata))

  # truncate times beyond tau to tau
  mydata$time[mydata$time > tau] <- tau

  mydata
}


#' Plot the estimated cumulative baseline mean function
#'
#' @param x A \code{wnpmle} object.
#' @param conf_int Logical; plot pointwise 95\% confidence bands (default
#'   \code{TRUE} if SE is available).
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Additional arguments passed to \code{plot}.
#' @export
plot.wnpmle <- function(x, conf_int = !anyNA(x$se_Lambda),
                         xlab = "Time",
                         ylab = expression(hat(Lambda)(t)),
                         main = "Estimated cumulative baseline mean",
                         ...) {
  t_obs  <- x$event_times
  Lambda <- x$Lambda

  plot(t_obs, Lambda, type = "s",
       xlab = xlab, ylab = ylab, main = main, ...)

  if (conf_int && !anyNA(x$se_Lambda)) {
    lwr <- pmax(Lambda - 1.96 * x$se_Lambda, 0)
    upr <- Lambda + 1.96 * x$se_Lambda
    lines(t_obs, lwr, type = "s", lty = 2, col = "grey50")
    lines(t_obs, upr, type = "s", lty = 2, col = "grey50")
  }
  invisible(x)
}
