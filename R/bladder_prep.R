#' Prepare bladder cancer data for wnpmle analysis
#'
#' Converts the \code{bladder1} dataset from the \pkg{survival} package into
#' the long format required by \code{\link{wnpmle_fit}}, with status codes
#' 0 (censored), 1 (recurrence), and 2 (death/dropout as terminal event).
#'
#' @param max_time Optional truncation time (default: maximum observed time).
#' @return A data frame in long format with columns \code{id}, \code{time},
#'   \code{status}, \code{rx}, \code{size}, and \code{number}.
#'
#' @details
#' The bladder cancer trial (Veterans Administration Cooperative Urological
#' Research Group) randomised patients to placebo or thiotepa. Recurrences of
#' bladder tumours are the recurrent event; the last observation per patient
#' is treated as a terminal event (status 2) if the patient left the study
#' for reasons other than a recurrence, otherwise as censored (status 0).
#'
#' @references
#' Byar, D.P. (1980). The Veterans Administration study of chemoprophylaxis
#' for recurrent stage I bladder tumors: comparisons of placebo, pyridoxine,
#' and thiotepa. In \emph{Bladder Tumors and Other Topics in Urological
#' Oncology}, 363-370. Plenum, New York.
#'
#' @examples
#' bdata <- bladder_prep()
#' head(bdata)
#' table(bdata$status)
#'
#' @export
bladder_prep <- function(max_time = NULL) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for bladder_prep().")

  blad <- survival::bladder1
  blad <- blad[order(blad$id, blad$stop), ]

  # build long format
  ids <- unique(blad$id)
  out <- lapply(ids, function(i) {
    sub <- blad[blad$id == i, , drop = FALSE]

    rows <- list()

    # recurrent events: event == 1, not the last observation
    for (j in seq_len(nrow(sub) - 1)) {
      if (sub$event[j] == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id     = i,
          time   = sub$stop[j],
          status = 1L,
          rx     = sub$rx[j],
          size   = sub$size[j],
          number = sub$number[j],
          stringsAsFactors = FALSE
        )
      }
    }

    # last row determines terminal/censored status
    ter <- sub[nrow(sub), , drop = FALSE]
    last_event <- ter$event[1]
    # event == 1: last recurrence (treat as recurrent + censored terminal)
    # event == 2: death/removal (terminal)
    # event == 0: censored
    if (last_event == 1L) {
      # last event is a recurrence — add it as recurrent, then add censored row
      rows[[length(rows) + 1]] <- data.frame(
        id     = i,
        time   = ter$stop[1],
        status = 1L,
        rx     = ter$rx[1],
        size   = ter$size[1],
        number = ter$number[1],
        stringsAsFactors = FALSE
      )
      term_status <- 0L
    } else {
      term_status <- if (last_event == 2L) 2L else 0L
    }

    rows[[length(rows) + 1]] <- data.frame(
      id     = i,
      time   = ter$stop[1],
      status = term_status,
      rx     = ter$rx[1],
      size   = ter$size[1],
      number = ter$number[1],
      stringsAsFactors = FALSE
    )

    do.call(rbind, rows)
  })

  out <- do.call(rbind, out)
  out$id <- as.integer(factor(out$id))

  if (!is.null(max_time))
    out <- out[out$time <= max_time, ]

  out[order(out$time), ]
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
