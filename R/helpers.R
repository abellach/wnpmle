# Internal helpers — not exported
# These are called by wnpmle_fit() for variance estimation.
# Score contributions use subject-level indexing (n x pscore),
# consistent with the simulation study implementation.

# ---- Score contributions (gradi matrix) ----
# Returns n x (numcov + num1) matrix indexed by subject id.
.compute_score <- function(model, rho, numcov, num1, n, n02, num2,
                            cov1, cov2, cov02, covc, beta, lambda, Lambda,
                            Lamc, Lam2, wnew, M1, M2, Mc) {

  gradi <- matrix(0, nrow = n, ncol = numcov + num1)

  ez1 <- as.numeric(exp(cov1  %*% beta))
  ez2 <- as.numeric(exp(cov2  %*% beta))
  ezc <- as.numeric(exp(covc  %*% beta))

  if (model == "boxcox") {

    # -- censored subjects --
    MGcbet <- (1 + ezc * Lamc)^(rho - 1)
    if (nrow(Mc) > 0) {
      for (j in seq_len(nrow(Mc))) {
        idj <- Mc$id[j]
        for (k in seq_len(numcov))
          gradi[idj, k] <- gradi[idj, k] +
            covc[j, k] * ezc[j] * Lamc[j] * MGcbet[j]
        for (k in seq_len(num1))
          gradi[idj, numcov + k] <- gradi[idj, numcov + k] +
            (M1$ind[k] <= Mc$ind[j]) * ezc[j] * MGcbet[j]
      }
    }

    # -- terminal events --
    if (num2 > 0) {
      Mglami   <- t(matrix(rep(lambda, num2), nrow = num1))
      MGLami   <- t(matrix(rep(Lambda, num2), nrow = num1))
      MGgrad2  <- (1 + ez2 * Lam2)^(rho - 1)
      vgrbet   <- ez2
      Mgcov2b  <- matrix(rep(vgrbet, num1), nrow = num2)
      Mgradin  <- Mgcov2b * MGLami
      MG1grad2 <- (1 + Mgradin)^(rho - 1)
      MG2grad2 <- (rho - 1) * (1 + Mgradin)^(rho - 2)
      MGgra1   <- wnew * Mglami * (MG1grad2 + MGLami * Mgcov2b * MG2grad2)
      Mgrb1    <- wnew * MG1grad2
      Mgrb2    <- wnew * MG2grad2 * Mglami

      for (j in seq_len(num2)) {
        idj <- M2$id[j]
        for (l in seq_len(numcov))
          gradi[idj, l] <- gradi[idj, l] +
            cov2[j, l] * ez2[j] * Lam2[j] * MGgrad2[j] +
            cov2[j, l] * vgrbet[j] * sum(MGgra1[j, ])
        if (num1 > 1) {
          for (l in seq_len(num1 - 1))
            gradi[idj, numcov + l] <- gradi[idj, numcov + l] +
              (M1$ind[l] <= M2$ind[j]) * ez2[j] * MGgrad2[j] +
              vgrbet[j] * Mgrb1[j, l] +
              (vgrbet[j]^2) * sum(Mgrb2[j, l:num1, drop = FALSE])
        }
        gradi[idj, numcov + num1] <- gradi[idj, numcov + num1] +
          (M1$ind[num1] <= M2$ind[j]) * ez2[j] * MGgrad2[j] +
          vgrbet[j] * Mgrb1[j, num1] +
          (vgrbet[j]^2) * Mgrb2[j, num1]
      }
    }

    # -- recurrent events --
    MGgrad1 <- (rho - 1) / (1 + ez1 * Lambda)
    for (i in seq_len(num1)) {
      idi <- M1$id[i]
      for (k in seq_len(numcov))
        gradi[idi, k] <- gradi[idi, k] -
          cov1[i, k] * (1 + ez1[i] * Lambda[i] * MGgrad1[i])
      for (k in seq_len(num1))
        gradi[idi, numcov + k] <- gradi[idi, numcov + k] -
          ((i == k) / lambda[k] + (k <= i) * ez1[i] * MGgrad1[i])
    }

  } else {  # logarithmic: G'(x) = 1/(1+rx),  G''(x) = -r/(1+rx)^2

    # -- censored subjects --
    MGcbet <- 1 / (1 + rho * ezc * Lamc)
    if (nrow(Mc) > 0) {
      for (j in seq_len(nrow(Mc))) {
        idj <- Mc$id[j]
        for (k in seq_len(numcov))
          gradi[idj, k] <- gradi[idj, k] +
            covc[j, k] * ezc[j] * Lamc[j] * MGcbet[j]
        for (k in seq_len(num1))
          gradi[idj, numcov + k] <- gradi[idj, numcov + k] +
            (M1$ind[k] <= Mc$ind[j]) * ezc[j] * MGcbet[j]
      }
    }

    # -- terminal events --
    if (num2 > 0) {
      Mglami   <- t(matrix(rep(lambda, num2), nrow = num1))
      MGLami   <- t(matrix(rep(Lambda, num2), nrow = num1))
      MGgrad2  <- 1 / (1 + rho * ez2 * Lam2)
      vgrbet   <- ez2
      Mgcov2b  <- matrix(rep(vgrbet, num1), nrow = num2)
      Mgradin  <- Mgcov2b * MGLami
      MG1grad2 <- 1    / (1 + rho * Mgradin)
      MG2grad2 <- (-rho) / (1 + rho * Mgradin)^2
      MGgra1   <- wnew * Mglami * (MG1grad2 + MGLami * Mgcov2b * MG2grad2)
      Mgrb1    <- wnew * MG1grad2
      Mgrb2    <- wnew * MG2grad2 * Mglami

      for (j in seq_len(num2)) {
        idj <- M2$id[j]
        for (l in seq_len(numcov))
          gradi[idj, l] <- gradi[idj, l] +
            cov2[j, l] * ez2[j] * Lam2[j] * MGgrad2[j] +
            cov2[j, l] * vgrbet[j] * sum(MGgra1[j, ])
        if (num1 > 1) {
          for (l in seq_len(num1 - 1))
            gradi[idj, numcov + l] <- gradi[idj, numcov + l] +
              (M1$ind[l] <= M2$ind[j]) * ez2[j] * MGgrad2[j] +
              vgrbet[j] * Mgrb1[j, l] +
              (vgrbet[j]^2) * sum(Mgrb2[j, l:num1, drop = FALSE])
        }
        gradi[idj, numcov + num1] <- gradi[idj, numcov + num1] +
          (M1$ind[num1] <= M2$ind[j]) * ez2[j] * MGgrad2[j] +
          vgrbet[j] * Mgrb1[j, num1] +
          (vgrbet[j]^2) * Mgrb2[j, num1]
      }
    }

    # -- recurrent events --
    MGgrad1 <- (-rho) / (1 + rho * ez1 * Lambda)
    for (i in seq_len(num1)) {
      idi <- M1$id[i]
      for (k in seq_len(numcov))
        gradi[idi, k] <- gradi[idi, k] -
          cov1[i, k] * (1 + ez1[i] * Lambda[i] * MGgrad1[i])
      for (k in seq_len(num1))
        gradi[idi, numcov + k] <- gradi[idi, numcov + k] -
          ((i == k) / lambda[k] + (k <= i) * ez1[i] * MGgrad1[i])
    }
  }

  gradi
}


# ---- Censoring correction (psi_subj matrix) ----
# Returns n x (numcov + num1) matrix of censoring-correction terms,
# accumulated at the subject level (matching simulation study implementation).
.censoring_correction <- function(model, rho, numcov, num1, n, n02, num2,
                                   cov2, beta, lambda, Lambda, wnew,
                                   M1, M2, M02) {

  ez2 <- as.numeric(exp(cov2 %*% beta))

  MLam    <- t(matrix(rep(Lambda, num2), nrow = num1))
  Mlam    <- t(matrix(rep(lambda, num2), nrow = num1))
  Mlamn02 <- matrix(rep(lambda, n02), nrow = num1)   # num1 x n02, row k = lambda[k]

  Mzbetn02 <- matrix(rep(ez2, n02), nrow = num2)     # num2 x n02, row j = ez2[j]
  Mzbet    <- matrix(rep(ez2, num1),   nrow = num2)
  MzbL     <- Mzbet * MLam

  if (model == "boxcox") {
    MG1 <- (1 + MzbL)^(rho - 1)
    MG2 <- (rho - 1) * (1 + MzbL)^(rho - 2)
  } else {
    MG1 <- 1    / (1 + rho * MzbL)
    MG2 <- (-rho) / (1 + rho * MzbL)^2
  }

  # Nelson-Aalen-type censoring hazard increments
  Mus1 <- matrix(0, n02, n02); diag(Mus1) <- M02$status0
  hazn <- (M02$status == 0) / (n02 - M02$idM02 + 1)
  Mus2 <- matrix(0, n02, n02)
  for (i in seq_len(n02)) Mus2[i, ] <- hazn * (M02$idM02[i] >= M02$idM02)
  Mu <- Mus1 - Mus2

  Mjk     <- wnew * (MG1 + Mzbet * MLam * MG2) * Mlam
  Muk.ind <- matrix(0, n02, num1)
  for (u in seq_len(n02)) Muk.ind[u, ] <- (M02$ind[u] <= M1$ind)

  Mju <- matrix(0, num2, n02)
  for (j in seq_len(num2))
    for (u in seq_len(n02))
      Mju[j, u] <- sum(Mjk[j, ] * Muk.ind[u, ])

  Mju.ind <- matrix(0, num2, n02)
  for (u in seq_len(n02)) Mju.ind[, u] <- (M2$ind <= M02$ind[u])
  Mju.new <- Mju * Mju.ind * Mzbetn02

  q1 <- matrix(0, numcov, n02)
  for (l in seq_len(numcov))
    for (u in seq_len(n02))
      q1[l, u] <- sum(cov2[, l] * Mju.new[, u])

  M2jl     <- wnew * Mzbet * MG1
  M2lu.ind <- t(Muk.ind)
  M2lu     <- matrix(0, num1, n02)
  for (l in seq_len(num1))
    for (u in seq_len(n02))
      M2lu[l, u] <- sum(M2jl[, l] * Mju.ind[, u])
  q21 <- M2lu.ind * M2lu

  M3jk <- (Mzbet^2) * wnew * MG2
  M3ku <- matrix(0, num1, n02)
  for (k in seq_len(num1))
    for (u in seq_len(n02))
      M3ku[k, u] <- (M02$ind[u] <= M1$ind[k]) *
                    sum(M3jk[, k] * (M2$ind <= M02$ind[u]))
  M3ku.new <- M3ku * Mlamn02

  M3kl.ind <- matrix(0, num1, num1)
  for (l in seq_len(num1)) M3kl.ind[l:num1, l] <- 1

  q22 <- matrix(0, num1, n02)
  for (l in seq_len(num1))
    for (u in seq_len(n02))
      q22[l, u] <- sum(M3ku.new[, u] * M3kl.ind[, l])

  q  <- rbind(q1, q21 + q22)
  pi <- n02 - M02$idM02 + 1

  # psi at the M02-row level
  psi_row <- matrix(0, nrow = n02, ncol = numcov + num1)
  for (k in seq_len(numcov + num1))
    for (i in seq_len(n02))
      psi_row[i, k] <- sum((q[k, ] / pi) * Mu[i, ])

  # accumulate back to subjects (each subject may have multiple M02 rows)
  psi_subj <- matrix(0, nrow = n, ncol = numcov + num1)
  for (i in seq_len(n02))
    psi_subj[M02$id[i], ] <- psi_subj[M02$id[i], ] + psi_row[i, ]

  psi_subj
}

# Folder for the compiled TMB templates.
# - options(wnpmle.cache_dir = "...") or the environment variable
#   WNPMLE_CACHE_DIR: use that folder (e.g. on a cluster, so that all jobs
#   share one compiled copy).
# - interactive sessions: tools::R_user_dir("wnpmle", "cache"), so the
#   templates are compiled only once per installed version of wnpmle and TMB.
#   Folders from older versions are removed.
# - otherwise (e.g. R CMD check): a folder in tempdir(), removed with the session.
.wnpmle_tmb_dir <- function() {
  user_dir <- getOption("wnpmle.cache_dir", Sys.getenv("WNPMLE_CACHE_DIR", ""))
  version_tag <- paste0("wnpmle", utils::packageVersion("wnpmle"),
                        "_TMB", utils::packageVersion("TMB"))
  if (nzchar(user_dir)) {
    dir <- file.path(path.expand(user_dir), version_tag)
  } else if (interactive()) {
    base <- tools::R_user_dir("wnpmle", which = "cache")
    dir  <- file.path(base, version_tag)
    if (dir.exists(base)) {
      old <- setdiff(list.dirs(base, recursive = FALSE), dir)
      if (length(old)) unlink(old, recursive = TRUE)
    }
  } else {
    dir <- file.path(tempdir(), "wnpmle_tmb")
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

# Compile and load the TMB template needed for the requested model.
# The package ships the .cpp templates in inst/tmb and compiles them lazily
# (see .wnpmle_tmb_dir() for where), avoiding compilation during library(wnpmle).
# part = "fit" (likelihood), "score" (subject-wise scores) or
# "corr" (censoring correction).
.wnpmle_compile_tmb_if_needed <- function(model, silent = TRUE, part = "fit") {
  suffix <- switch(model,
    boxcox = "BC_tmb",
    log    = "log_tmb",
    stop("Unknown model: ", model, call. = FALSE)
  )
  prefix <- switch(part,
    fit   = "fn_",
    score = "score_",
    corr  = "corr_",
    stop("Unknown part: ", part, call. = FALSE)
  )
  dll_base <- paste0(prefix, suffix)
  cpp_name <- paste0(dll_base, ".cpp")
  cpp_src <- system.file("tmb", cpp_name, package = "wnpmle", mustWork = TRUE)

  build_dir <- .wnpmle_tmb_dir()

  cpp_tmp <- file.path(build_dir, cpp_name)
  dll_tmp <- TMB::dynlib(file.path(build_dir, dll_base))

  if (!file.exists(dll_tmp) || file.info(dll_tmp)$mtime < file.info(cpp_src)$mtime) {
    file.copy(cpp_src, cpp_tmp, overwrite = TRUE)
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(build_dir)
    TMB::compile(cpp_name, silent = silent)
  }

  if (!dll_base %in% names(getLoadedDLLs())) {
    dyn.load(dll_tmp)
  }

  dll_base
}


# ---- Score contributions via TMB (score_BC_tmb.cpp / score_log_tmb.cpp) ----
# Returns the numi x (numcov + num1) matrix of subject-wise derivatives of the
# negative log-likelihood with respect to (beta, lambda); identical to
# .compute_score() (checked to 1e-13).
.compute_score_tmb <- function(model, data_tmb, parameters, par_hat,
                               numcov, num1, numi, lambda, M1, M2, M02,
                               silent = TRUE) {
  dll <- .wnpmle_compile_tmb_if_needed(model, silent = silent, part = "score")
  d <- data_tmb
  d$id1  <- as.integer(M1$id  - 1L)
  d$id02 <- as.integer(M02$id - 1L)
  d$id2  <- as.integer(M2$id  - 1L)
  d$numi <- as.integer(numi)
  obj <- TMB::MakeADFun(d, parameters, DLL = dll, ADreport = TRUE, silent = silent)
  J <- obj$gr(par_hat)                                  # wrt (beta, alpha)
  idx <- numcov + seq_len(num1)
  J[, idx] <- sweep(J[, idx, drop = FALSE], 2, lambda, "/")   # alpha -> lambda
  J
}

# ---- Censoring correction via TMB (corr_BC_tmb.cpp / corr_log_tmb.cpp) ----
# Returns psi_subj (numi x (numcov + num1)) with the same sign as
# .censoring_correction(); it enters the meat as  gradi - psi_subj.
.censoring_correction_tmb <- function(model, data_tmb, parameters, par_hat,
                                      numi, M1, M2, M02, silent = TRUE) {
  dll <- .wnpmle_compile_tmb_if_needed(model, silent = silent, part = "corr")
  d <- data_tmb
  d$ind1     <- as.integer(M1$ind)
  d$ind2     <- as.integer(M2$ind)
  d$ind02    <- as.integer(M02$ind)
  d$status02 <- as.integer(M02$status0)
  d$id02     <- as.integer(M02$id - 1L)
  d$numi     <- as.integer(numi)
  obj <- TMB::MakeADFun(d, parameters, DLL = dll, silent = silent)
  obj$report(par_hat)$psi_subj
}
