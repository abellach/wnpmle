// corr_BC_tmb.cpp
//
// Censoring correction psi_subj (kappa_i in the paper) for the Box-Cox
// recurrent-event model. Line-by-line translation of the R code in
// sim_rec_BCrho2tau5_n100_10k_cluster.R, lines 324-386 (same names:
// Mjk, Mju, q1, q21, q22, q, pi, hazn, Mu, psi_row, psi_subj).
//
// Returns psi_subj with the SAME sign as the R code, so the two can be
// compared directly. In the meat use  gradi - psi_subj  (see check of 7 Oct).
//
// Usage in R (at the fitted parameters):
//   obj_c <- MakeADFun(data_corr, parameters, DLL = "corr_BC_tmb", silent = TRUE)
//   psi_subj <- obj_c$report(opt$par)$psi_subj
//
// The computation only runs in report(); the objective itself is 0.
//
// Data: the same as fn_BC_tmb.cpp, plus
//   ind1  (num1): M1$ind,  row position of each recurrent event
//   ind2  (num2): M2$ind,  row position of each terminal event
//   ind02 (n02):  M02$ind, row position of each censored/terminal row
//   status02 (n02): M02$status0 (1 = censored, 0 = terminal)
//   id02  (n02):  subject of each censored/terminal row, 0-based (M02$id - 1)
//   numi:         number of subjects
//
#define TMB_LIB_INIT R_init_corr_BC_tmb
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(cov1);
  DATA_MATRIX(cov2);
  DATA_MATRIX(cov02);
  DATA_IVECTOR(idx02);
  DATA_IVECTOR(idx2);
  DATA_VECTOR(kmc1);
  DATA_VECTOR(kmc2);
  DATA_SCALAR(rho);
  DATA_IVECTOR(ind1);
  DATA_IVECTOR(ind2);
  DATA_IVECTOR(ind02);
  DATA_IVECTOR(status02);
  DATA_IVECTOR(id02);
  DATA_INTEGER(numi);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

  if (isDouble<Type>::value) {

    const int num1   = cov1.rows();
    const int num2   = cov2.rows();
    const int n02    = cov02.rows();
    const int numcov = cov1.cols();
    const int pscore = numcov + num1;

    vector<Type> lambda = exp(alpha);
    vector<Type> Lambda(num1);
    Type csum = Type(0);
    for (int k = 0; k < num1; ++k) { csum += lambda(k); Lambda(k) = csum; }

    vector<Type> eta2 = cov2 * beta;
    vector<Type> e2   = exp(eta2);

    // ---- wnew, MG1, MG2 (num2 x num1) ----
    // wnew[j,k] = (M1$ind[k] > M2$ind[j]) * M1$kmc[k] / M2$kmc[j]
    matrix<Type> wnew(num2, num1), MG1(num2, num1), MG2(num2, num1);
    for (int j = 0; j < num2; ++j) {
      for (int k = 0; k < num1; ++k) {
        wnew(j, k) = (ind1(k) > ind2(j)) ? kmc1(k) / kmc2(j) : Type(0);
        Type MzbL  = e2(j) * Lambda(k);
        MG1(j, k)  = pow(Type(1) + MzbL, rho - Type(1));
        MG2(j, k)  = (rho - Type(1)) * pow(Type(1) + MzbL, rho - Type(2));
      }
    }

    // ---- q (pscore x n02) = rbind(q1, q21 + q22) ----
    matrix<Type> q(pscore, n02);
    q.setZero();

    for (int u = 0; u < n02; ++u) {
      for (int j = 0; j < num2; ++j) {
        if (!(ind2(j) <= ind02(u))) continue;          // Mju.ind[j,u]

        // Mju[j,u] = sum_k Mjk[j,k] * Muk.ind[u,k],
        // Mjk = wnew * (MG1 + Mzbet * MLam * MG2) * Mlam
        Type Mju = Type(0);
        for (int k = 0; k < num1; ++k) {
          if (ind02(u) <= ind1(k)) {                   // Muk.ind[u,k]
            Mju += wnew(j, k) * (MG1(j, k) + e2(j) * Lambda(k) * MG2(j, k)) * lambda(k);
          }
        }
        // q1[l,u] = sum_j cov2[j,l] * Mju[j,u] * Mju.ind[j,u] * Mzbet[j]
        for (int l = 0; l < numcov; ++l) q(l, u) += cov2(j, l) * Mju * e2(j);

        // q21[l,u] = Muk.ind[u,l] * sum_j M2jl[j,l] * Mju.ind[j,u],
        // M2jl = wnew * Mzbet * MG1
        for (int l = 0; l < num1; ++l) {
          if (ind02(u) <= ind1(l)) {
            q(numcov + l, u) += wnew(j, l) * e2(j) * MG1(j, l);
          }
        }
      }

      // q22[l,u] = sum_{k >= l} M3ku.new[k,u],
      // M3ku[k,u] = Muk.ind[u,k] * sum_j M3jk[j,k] * Mju.ind[j,u],
      // M3jk = Mzbet^2 * wnew * MG2,  M3ku.new = M3ku * lambda[k]
      Type suffix = Type(0);
      for (int k = num1 - 1; k >= 0; --k) {
        Type M3ku = Type(0);
        if (ind02(u) <= ind1(k)) {
          for (int j = 0; j < num2; ++j) {
            if (ind2(j) <= ind02(u)) M3ku += e2(j) * e2(j) * wnew(j, k) * MG2(j, k);
          }
        }
        suffix += M3ku * lambda(k);
        q(numcov + k, u) += suffix;
      }
    }

    // ---- psi_row[i,k] = sum_u (q[k,u] / pi[u]) * Mu[i,u] ----
    // pi[u] = n02 - idM02[u] + 1,  hazn[u] = status0[u] / pi[u],
    // Mu[i,u] = status0[i] * (i == u) - hazn[u] * (u <= i)
    // (the sum over u <= i is accumulated as a running sum, same result)
    matrix<Type> psi_row(n02, pscore);
    vector<Type> running(pscore);
    running.setZero();
    for (int i = 0; i < n02; ++i) {
      Type pi_i   = Type(n02 - i);
      Type hazn_i = Type(status02(i)) / pi_i;
      for (int k = 0; k < pscore; ++k) {
        running(k) += hazn_i * q(k, i) / pi_i;
        psi_row(i, k) = Type(status02(i)) * q(k, i) / pi_i - running(k);
      }
    }

    // ---- psi_subj: add rows to subjects ----
    matrix<Type> psi_subj(numi, pscore);
    psi_subj.setZero();
    for (int i = 0; i < n02; ++i)
      for (int k = 0; k < pscore; ++k)
        psi_subj(id02(i), k) += psi_row(i, k);

    REPORT(q);
    REPORT(psi_subj);
  }

  return Type(0);
}
