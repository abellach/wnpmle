# Internal helpers — not exported
# These are called by wnpmle_fit() for variance estimation.
# Score contributions use subject-level indexing (n x pscore),
# consistent with the simulation study implementation.

# ======================================================================
# OLD VERSIONS (loop-based) — kept for validation only
# Will be removed once sandwich_adj and sandwich_adjold are confirmed equal
# ======================================================================

.compute_score_old <- function(model, rho, numcov, num1, n, n02, num2,
                                cov1, cov2, cov02, covc, beta, lambda, Lambda,
                                Lamc, Lam2, wnew, M1, M2, Mc) {

  gradi <- matrix(0, nrow = n, ncol = numcov + num1)

  ez1 <- as.numeric(exp(cov1  %*% beta))
  ez2 <- as.numeric(exp(cov2  %*% beta))
  ezc <- as.numeric(exp(covc  %*% beta))

  if (model == "boxcox") {

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

  } else {

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


.censoring_correction_old <- function(model, rho, numcov, num1, n, n02, num2,
                                       cov2, beta, lambda, Lambda, wnew,
                                       M1, M2, M02) {

  ez2 <- as.numeric(exp(cov2 %*% beta))

  MLam    <- t(matrix(rep(Lambda, num2), nrow = num1))
  Mlam    <- t(matrix(rep(lambda, num2), nrow = num1))
  Mlamn02 <- t(matrix(rep(lambda, n02),  nrow = n02))
  Mzbetn02 <- t(matrix(rep(ez2, n02),  nrow = n02))
  Mzbet    <- matrix(rep(ez2, num1),   nrow = num2)
  MzbL     <- Mzbet * MLam

  if (model == "boxcox") {
    MG1 <- (1 + MzbL)^(rho - 1)
    MG2 <- (rho - 1) * (1 + MzbL)^(rho - 2)
  } else {
    MG1 <- 1    / (1 + rho * MzbL)
    MG2 <- (-rho) / (1 + rho * MzbL)^2
  }

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

  psi_row <- matrix(0, nrow = n02, ncol = numcov + num1)
  for (k in seq_len(numcov + num1))
    for (i in seq_len(n02))
      psi_row[i, k] <- sum((q[k, ] / pi) * Mu[i, ])

  psi_subj <- matrix(0, nrow = n, ncol = numcov + num1)
  for (i in seq_len(n02))
    psi_subj[M02$id[i], ] <- psi_subj[M02$id[i], ] + psi_row[i, ]

  psi_subj
}


# ======================================================================
# NEW VERSIONS (vectorised) — matrix operations replacing nested loops
# ======================================================================

.compute_score <- function(model, rho, numcov, num1, n, n02, num2,
                            cov1, cov2, cov02, covc, beta, lambda, Lambda,
                            Lamc, Lam2, wnew, M1, M2, Mc) {

  gradi <- matrix(0, nrow = n, ncol = numcov + num1)

  ez1 <- as.numeric(exp(cov1  %*% beta))
  ez2 <- as.numeric(exp(cov2  %*% beta))
  ezc <- as.numeric(exp(covc  %*% beta))

  if (model == "boxcox") {

    MGcbet <- (1 + ezc * Lamc)^(rho - 1)
    if (nrow(Mc) > 0) {
      gradi[Mc$id, seq_len(numcov)] <- gradi[Mc$id, seq_len(numcov)] +
        (ezc * Lamc * MGcbet) * covc
      ind_c1 <- outer(Mc$ind, M1$ind, "<=")
      gradi[Mc$id, numcov + seq_len(num1)] <-
        gradi[Mc$id, numcov + seq_len(num1)] +
        ind_c1 * (ezc * MGcbet)
    }

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

      beta_term1 <- (ez2 * Lam2 * MGgrad2) * cov2
      beta_term2 <- vgrbet * rowSums(MGgra1) * cov2
      gradi[M2$id, seq_len(numcov)] <-
        gradi[M2$id, seq_len(numcov)] + beta_term1 + beta_term2

      ind_21    <- outer(M1$ind, M2$ind, "<=")
      lam_term1 <- t(ind_21) * (ez2 * MGgrad2)
      lam_term2 <- vgrbet * Mgrb1
      Mgrb2_cumrev <- t(apply(Mgrb2, 1, function(x) rev(cumsum(rev(x)))))
      lam_term3 <- (vgrbet^2) * Mgrb2_cumrev
      gradi[M2$id, numcov + seq_len(num1)] <-
        gradi[M2$id, numcov + seq_len(num1)] +
        lam_term1 + lam_term2 + lam_term3
    }

    MGgrad1 <- (rho - 1) / (1 + ez1 * Lambda)
    gradi[M1$id, seq_len(numcov)] <- gradi[M1$id, seq_len(numcov)] -
      (1 + ez1 * Lambda * MGgrad1) * cov1
    diag_term  <- matrix(0, num1, num1); diag(diag_term) <- -1 / lambda
    ind_11     <- outer(seq_len(num1), seq_len(num1), ">=")
    lower_term <- -ind_11 * (ez1 * MGgrad1)
    gradi[M1$id, numcov + seq_len(num1)] <-
      gradi[M1$id, numcov + seq_len(num1)] + diag_term + lower_term

  } else {

    MGcbet <- 1 / (1 + rho * ezc * Lamc)
    if (nrow(Mc) > 0) {
      gradi[Mc$id, seq_len(numcov)] <- gradi[Mc$id, seq_len(numcov)] +
        (ezc * Lamc * MGcbet) * covc
      ind_c1 <- outer(Mc$ind, M1$ind, "<=")
      gradi[Mc$id, numcov + seq_len(num1)] <-
        gradi[Mc$id, numcov + seq_len(num1)] +
        ind_c1 * (ezc * MGcbet)
    }

    if (num2 > 0) {
      Mglami   <- t(matrix(rep(lambda, num2), nrow = num1))
      MGLami   <- t(matrix(rep(Lambda, num2), nrow = num1))
      MGgrad2  <- 1 / (1 + rho * ez2 * Lam2)
      vgrbet   <- ez2
      Mgcov2b  <- matrix(rep(vgrbet, num1), nrow = num2)
      Mgradin  <- Mgcov2b * MGLami
      MG1grad2 <- 1     / (1 + rho * Mgradin)
      MG2grad2 <- (-rho) / (1 + rho * Mgradin)^2
      MGgra1   <- wnew * Mglami * (MG1grad2 + MGLami * Mgcov2b * MG2grad2)
      Mgrb1    <- wnew * MG1grad2
      Mgrb2    <- wnew * MG2grad2 * Mglami

      beta_term1 <- (ez2 * Lam2 * MGgrad2) * cov2
      beta_term2 <- vgrbet * rowSums(MGgra1) * cov2
      gradi[M2$id, seq_len(numcov)] <-
        gradi[M2$id, seq_len(numcov)] + beta_term1 + beta_term2

      ind_21    <- outer(M1$ind, M2$ind, "<=")
      lam_term1 <- t(ind_21) * (ez2 * MGgrad2)
      lam_term2 <- vgrbet * Mgrb1
      Mgrb2_cumrev <- t(apply(Mgrb2, 1, function(x) rev(cumsum(rev(x)))))
      lam_term3 <- (vgrbet^2) * Mgrb2_cumrev
      gradi[M2$id, numcov + seq_len(num1)] <-
        gradi[M2$id, numcov + seq_len(num1)] +
        lam_term1 + lam_term2 + lam_term3
    }

    MGgrad1 <- (-rho) / (1 + rho * ez1 * Lambda)
    gradi[M1$id, seq_len(numcov)] <- gradi[M1$id, seq_len(numcov)] -
      (1 + ez1 * Lambda * MGgrad1) * cov1
    diag_term  <- matrix(0, num1, num1); diag(diag_term) <- -1 / lambda
    ind_11     <- outer(seq_len(num1), seq_len(num1), ">=")
    lower_term <- -ind_11 * (ez1 * MGgrad1)
    gradi[M1$id, numcov + seq_len(num1)] <-
      gradi[M1$id, numcov + seq_len(num1)] + diag_term + lower_term
  }

  gradi
}


.censoring_correction <- function(model, rho, numcov, num1, n, n02, num2,
                                   cov2, beta, lambda, Lambda, wnew,
                                   M1, M2, M02) {

  ez2 <- as.numeric(exp(cov2 %*% beta))

  MLam    <- t(matrix(rep(Lambda, num2), nrow = num1))
  Mlam    <- t(matrix(rep(lambda, num2), nrow = num1))
  Mlamn02 <- t(matrix(rep(lambda, n02),  nrow = n02))
  Mzbetn02 <- t(matrix(rep(ez2, n02),  nrow = n02))
  Mzbet    <- matrix(rep(ez2, num1),   nrow = num2)
  MzbL     <- Mzbet * MLam

  if (model == "boxcox") {
    MG1 <- (1 + MzbL)^(rho - 1)
    MG2 <- (rho - 1) * (1 + MzbL)^(rho - 2)
  } else {
    MG1 <- 1     / (1 + rho * MzbL)
    MG2 <- (-rho) / (1 + rho * MzbL)^2
  }

  Mus1 <- matrix(0, n02, n02); diag(Mus1) <- M02$status0
  hazn <- (M02$status == 0) / (n02 - M02$idM02 + 1)
  Mus2 <- outer(M02$idM02, M02$idM02, ">=") * hazn
  Mu   <- Mus1 - Mus2

  Muk.ind <- outer(M02$ind, M1$ind, "<=")
  Mju.ind <- outer(M2$ind,  M02$ind, "<=")

  Mjk     <- wnew * (MG1 + Mzbet * MLam * MG2) * Mlam
  Mju     <- Mjk %*% t(Muk.ind)
  Mju.new <- Mju * Mju.ind * Mzbetn02
  q1      <- t(cov2) %*% Mju.new

  M2jl <- wnew * Mzbet * MG1
  M2lu <- t(M2jl) %*% Mju.ind
  q21  <- t(Muk.ind) * M2lu

  M3jk     <- (Mzbet^2) * wnew * MG2
  M3ku     <- t(Muk.ind) * (t(M3jk) %*% Mju.ind)
  M3ku.new <- M3ku * Mlamn02
  q22      <- apply(M3ku.new, 2, function(col) rev(cumsum(rev(col))))

  q  <- rbind(q1, q21 + q22)
  pi <- n02 - M02$idM02 + 1

  psi_row  <- Mu %*% t(q / pi)

  psi_subj <- matrix(0, nrow = n, ncol = numcov + num1)
  for (i in seq_len(n02))
    psi_subj[M02$id[i], ] <- psi_subj[M02$id[i], ] + psi_row[i, ]

  psi_subj
}
