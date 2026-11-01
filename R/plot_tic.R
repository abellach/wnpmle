#' Takeuchi information criterion over the transformation parameter
#'
#' Fits the Box-Cox and logarithmic transformation models over grids of the
#' transformation parameter and plots Takeuchi's information criterion
#' (\code{\link{TIC}}) for both families on a single axis, in the same layout
#' as \code{\link{plot_loglik}}: the logarithmic parameter \eqn{r} on the left,
#' the Box-Cox parameter \eqn{\rho} on the right. The two families meet at
#' zero: \eqn{r \to 0} and \eqn{\rho = 1} give the Ghosh-Lin model, and
#' \eqn{\rho \to 0} and \eqn{r = 1} the proportional odds model.
#'
#' Every grid point needs a fit with the censoring-corrected sandwich and a
#' fit with the inverse Fisher information, so this is much slower than
#' \code{\link{plot_loglik}}; the default grids are therefore coarser.
#'
#' @inheritParams plot_loglik
#' @param rho_grid Box-Cox parameter values (default \code{seq(0.05, 1.2, by = 0.05)}).
#' @param r_grid Logarithmic parameter values (default \code{seq(0.05, 1.2, by = 0.05)}).
#' @param show_aic Logical; if \code{TRUE} (default), AIC is added in grey
#'   for comparison.
#'
#' @return A data frame with columns \code{model}, \code{param}, \code{loglik},
#'   \code{penalty}, \code{p}, \code{AIC} and \code{TIC}, invisibly.
#'
#' @seealso \code{\link{TIC}}, \code{\link{plot_loglik}}
#'
#' @examples
#' \dontrun{
#'   dat <- bladder_prep()
#'   res <- plot_tic(Surv(time, status) ~ treat + num + size, data = dat,
#'                   id = "id", rho_grid = seq(0.2, 1.2, by = 0.2),
#'                   r_grid = seq(0.2, 1.2, by = 0.2))
#'   res
#' }
#'
#' @export
plot_tic <- function(formula, data, id = "id",
                     rho_grid    = seq(0.05, 1.2, by = 0.05),
                     r_grid      = seq(0.05, 1.2, by = 0.05),
                     tau         = NULL,
                     mark_points = TRUE,
                     show_aic    = TRUE,
                     file        = NULL,
                     width       = 6,
                     height      = 4.5,
                     verbose     = TRUE,
                     ...) {

  # TIC (with details) for one model; NA row if a fit fails
  one_tic <- function(model, p_val, init_beta) {
    fs <- tryCatch(wnpmle_fit(formula, data = data, id = id, model = model,
                              rho = p_val, tau = tau, se = "sandwich_adj",
                              init_beta = init_beta, ...),
                   error = function(e) NULL)
    if (is.null(fs)) return(list(res = rep(NA_real_, 5), beta = init_beta))
    ff <- tryCatch(wnpmle_fit(formula, data = data, id = id, model = model,
                              rho = p_val, tau = tau, se = "fisher",
                              init_beta = fs$coefficients, ...),
                   error = function(e) NULL)
    if (is.null(ff)) return(list(res = rep(NA_real_, 5), beta = fs$coefficients))
    d <- tryCatch(suppressWarnings(TIC(fs, fisher = ff, details = TRUE)),
                  error = function(e) rep(NA_real_, 5))
    list(res = as.numeric(d[c("loglik", "penalty", "p", "AIC", "TIC")]),
         beta = fs$coefficients)
  }

  run_grid <- function(model, grid, label) {
    if (verbose) cat("Fitting", label, "grid (", length(grid), "models)...\n")
    out <- matrix(NA_real_, length(grid), 5,
                  dimnames = list(NULL, c("loglik", "penalty", "p", "AIC", "TIC")))
    init_beta <- NULL
    for (k in seq_along(grid)) {
      o <- one_tic(model, grid[k], init_beta)
      out[k, ]  <- o$res
      init_beta <- o$beta                     # warm start
      if (verbose && k %% 5 == 0) cat("  ", label, ":", k, "/", length(grid), "\n")
    }
    out
  }

  res_BC  <- run_grid("boxcox", rho_grid, "Box-Cox")
  res_log <- run_grid("log",    r_grid,   "Log")

  # join the families at 0 as in plot_loglik():
  # Box-Cox rho -> 0 is log r = 1, log r -> 0 is Box-Cox rho = 1
  i_r1_val   <- which.min(abs(r_grid   - 1))
  i_rho1_val <- which.min(abs(rho_grid - 1))
  tic_BC.new   <- c(res_log[i_r1_val, "TIC"],  res_BC[,  "TIC"])
  tic_log.new  <- c(res_BC[i_rho1_val, "TIC"], res_log[, "TIC"])
  aic_BC.new   <- c(res_log[i_r1_val, "AIC"],  res_BC[,  "AIC"])
  aic_log.new  <- c(res_BC[i_rho1_val, "AIC"], res_log[, "AIC"])
  rho_grid.new <- c(0, rho_grid)
  r_grid.new   <- c(0, r_grid)

  # ---- plot ----
  if (!is.null(file)) pdf(file, width = width, height = height, useDingbats = FALSE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  yvals    <- c(tic_log.new, tic_BC.new, if (show_aic) c(aic_log.new, aic_BC.new))
  ylim_all <- range(yvals, finite = TRUE)
  max_r    <- ceiling(max(r_grid.new)   / 0.4) * 0.4
  max_rho  <- ceiling(max(rho_grid.new) / 0.4) * 0.4
  xlim_all <- c(-max_r, max_rho)

  par(mar = c(5, 4.5, 3, 3), mgp = c(2.5, 0.5, 0), tcl = -0.25)
  plot(NA, xlim = xlim_all, ylim = ylim_all, xlab = "",
       ylab = if (show_aic) "TIC (black), AIC (grey)" else "TIC", axes = FALSE)
  box()

  if (show_aic) {
    lines(-r_grid.new,  aic_log.new, lwd = 1.5, lty = 2, col = "grey60")
    lines(rho_grid.new, aic_BC.new,  lwd = 1.5, col = "grey60")
  }
  lines(-r_grid.new,  tic_log.new, lwd = 2, lty = 2)
  lines(rho_grid.new, tic_BC.new,  lwd = 2)
  abline(v = 0, lty = 3, col = "grey60")

  ticks_left  <- round(seq(-max_r, 0,       by = 0.4), 2)
  ticks_right <- round(seq(0,      max_rho, by = 0.4), 2)
  axis(1, at = c(ticks_left, ticks_right[-1]),
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
    points(-r_grid.new[i_r1], tic_log.new[i_r1], pch = 16, cex = 1)
    text(-r_grid.new[i_r1], tic_log.new[i_r1], labels = "r = 1",
         pos = 4, offset = 0.7, xpd = NA)
    points(rho_grid.new[i_rho1], tic_BC.new[i_rho1], pch = 1, cex = 1.5)
    text(rho_grid.new[i_rho1], tic_BC.new[i_rho1], labels = expression(rho == 1),
         pos = 1, offset = 0.9, xpd = NA)
  }

  if (!is.null(file)) {
    dev.off()
    cat("Plot saved to", file, "\n")
  }

  # ---- report optima (smallest TIC) ----
  out <- data.frame(
    model = c(rep("log", length(r_grid)), rep("boxcox", length(rho_grid))),
    param = c(r_grid, rho_grid),
    rbind(res_log, res_BC),
    row.names = NULL
  )
  best <- function(m) {
    s <- out[out$model == m & is.finite(out$TIC), ]
    s[which.min(s$TIC), ]
  }
  b_bc <- best("boxcox"); b_lg <- best("log")
  cat("\nTransformation parameters with the smallest TIC:\n")
  print(data.frame(
    Model     = c("Box-Cox (Ghosh-Lin at rho=1)", "Log (Prop. odds at r=1)"),
    Parameter = c(paste("rho =", b_bc$param), paste("r =", b_lg$param)),
    TIC       = round(c(b_bc$TIC, b_lg$TIC), 3),
    AIC       = round(c(b_bc$AIC, b_lg$AIC), 3)
  ), row.names = FALSE)

  invisible(out)
}
