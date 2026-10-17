// corr_cmprsk_BC.hpp
//
// Correction of the sandwich for the estimated IPC weights (KM of the
// censoring distribution), Box-Cox competing-risks model. Same likelihood and
// data layout as fn_cmprsk_BC_tmb.cpp / score_cmprsk_BC_tmb.cpp.
//
// Technical report, Section B.4:
//   kappa_i = int q(u) / pi(u) dM_i^c(u)
//   Sigma   = sum_i (eta_i + kappa_i)^2
// Everything here is computed from the NEGATIVE log-likelihood (like gradi),
// and the linearization of KM-hat/KM carries a minus sign, so the corrected
// meat is sum_i (gradi_i - psi_i)^2 = sum_i (eta_i + kappa_i)^2.
// (The old code used gradi + psi; sign fixed 8-9 Oct.)
//
// Derivation used here (positions 0..n-1 in time order):
//   w(j,k) = KM(X_k)/KM(X_j) = prod over censorings u with p2_j < u < p1_k of
//            (1 - 1/Y_u),  Y_u = n - u  (number at risk at position u)
//   => d w(j,k) / d dLambda^c(u) = - w(j,k) D_jk(u),  D_jk(u) = 1(p2_j < u < p1_k)
//   The weights enter only through
//      f5 = sum_j sum_k w(j,k) h_jk,  h_jk = e_j lambda_k G'(e_j Lambda_k).
//   q(u) = sum_j sum_k D_jk(u) w(j,k) d h_jk / d theta,  theta = (beta, lambda):
//      d h_jk / d beta_l   = z_jl e_j lambda_k [G'(x_jk) + e_j Lambda_k G''(x_jk)]
//      d h_jk / d lambda_m = 1(m = k) e_j G'(x_jk) + 1(m <= k) e_j^2 lambda_k G''(x_jk)
//   psi_i = sum_u q(u) / Y_u * dM_i^c(u),
//      dM_i^c(u) = 1(i = u) - 1(i >= u) / Y_u     (u a censoring position)
//
// For fixed u, D_jk(u) = 1 exactly for j with p2_j < u and k >= k0(u), the
// first type-1 event after u. Processing u in increasing order and adding the
// type-2 subjects j as they enter keeps the cost at O(num2*num1*numcov +
// numc*num1) instead of the triple loop of the old code.
//
// Use: MakeADFun(..., silent = TRUE)$report() returns psi (n x P, columns
// (beta, lambda)) and q (P x numc). The objective itself is a dummy 0.
//

// G'(x)  = (1+x)^(rho-1)
//
// Model "corr_cmprsk_BC" of the wnpmle library (competing risks). Called from
// src/wnpmle.cpp when the data list passed to MakeADFun() contains
// model = "corr_cmprsk_BC". The code between the two marker lines is the stand-alone
// template verified on 16 Oct 2026 (test_step12/3/4.R); only its first line
// was changed (objective_function -> model).
#ifndef CORR_CMPRSK_BC_HPP
#define CORR_CMPRSK_BC_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

namespace wnpmle_corr_cmprsk_BC {
// ---- original template begins ----
template<class Type> Type G1f(Type x, Type rho) { return pow(Type(1) + x, rho - Type(1)); }
// G''(x) = (rho-1)(1+x)^(rho-2)
template<class Type> Type G2f(Type x, Type rho) { return (rho - Type(1)) * pow(Type(1) + x, rho - Type(2)); }

template<class Type>
Type model(objective_function<Type>* obj)
{
  DATA_MATRIX(cova);      // n x numcov, all subjects ordered by time
  DATA_IVECTOR(ind1);     // num1, 0-based positions of type-1 events (increasing)
  DATA_IVECTOR(ind2);     // num2, 0-based positions of type-2 events (increasing)
  DATA_IVECTOR(indc);     // numc, 0-based positions of censorings (increasing)
  DATA_MATRIX(w);         // num2 x num1, IPCW weights
  DATA_SCALAR(rho);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

  if (isDouble<Type>::value) {
    int n = cova.rows(), numcov = cova.cols();
    int num1 = ind1.size(), num2 = ind2.size(), numc = indc.size();
    int P = numcov + num1;

    vector<Type> lambda(num1), Lambda(num1);
    Type csum = Type(0);
    for (int k = 0; k < num1; k++) { lambda(k) = exp(alpha(k)); csum += lambda(k); Lambda(k) = csum; }

    // per type-2 subject j: a_jk (point term), and suffix sums over k of
    //   s1_jk = w e lambda_k (G' + e Lambda_k G'')  (for beta rows, times z_j)
    //   t_jk  = w lambda_k G''                       (for lambda rows, times e^2)
    matrix<Type> A(num2, num1), S1(num2, num1 + 1), T(num2, num1 + 1);
    vector<Type> e2(num2);
    for (int j = 0; j < num2; j++) {
      Type eta = Type(0);
      for (int l = 0; l < numcov; l++) eta += cova(ind2(j), l) * beta(l);
      Type e = exp(eta);
      e2(j) = e;
      S1(j, num1) = Type(0); T(j, num1) = Type(0);
      for (int k = num1 - 1; k >= 0; k--) {
        Type wjk = w(j, k);
        Type a = Type(0), s = Type(0), t = Type(0);
        if (wjk != Type(0)) {
          Type x  = e * Lambda(k);
          Type g1 = G1f(x, rho), g2 = G2f(x, rho);
          a = wjk * e * g1;
          s = wjk * e * lambda(k) * (g1 + e * Lambda(k) * g2);
          t = wjk * lambda(k) * g2;
        }
        A(j, k) = a;
        S1(j, k) = S1(j, k + 1) + s;
        T(j, k)  = T(j, k + 1) + t;
      }
    }

    // q(u) for each censoring u, accumulating type-2 subjects with p2_j < u
    matrix<Type> q(P, numc);
    q.setZero();
    vector<Type> accA(num1), accT(num1 + 1);
    matrix<Type> accS(numcov, num1 + 1);
    accA.setZero(); accT.setZero(); accS.setZero();
    int jnext = 0, k0 = 0;
    for (int c = 0; c < numc; c++) {
      int u = indc(c);
      while (jnext < num2 && ind2(jnext) < u) {
        int j = jnext;
        for (int k = 0; k < num1; k++) accA(k) += A(j, k);
        for (int k = 0; k <= num1; k++) {
          accT(k) += e2(j) * e2(j) * T(j, k);
          for (int l = 0; l < numcov; l++) accS(l, k) += cova(ind2(j), l) * S1(j, k);
        }
        jnext++;
      }
      while (k0 < num1 && ind1(k0) < u) k0++;   // first type-1 event after u
      if (k0 >= num1) continue;                  // no type-1 event after u: q(u) = 0
      for (int l = 0; l < numcov; l++) q(l, c) = accS(l, k0);
      for (int m = 0; m < num1; m++) {
        Type v = accT(m > k0 ? m : k0);
        if (m >= k0) v += accA(m);
        q(numcov + m, c) = v;
      }
    }

    // psi_i = sum_u q(u)/Y_u * (1(i = u) - 1(i >= u)/Y_u)
    matrix<Type> psi(n, P);
    psi.setZero();
    vector<Type> cum(P);
    cum.setZero();
    int c = 0;
    for (int i = 0; i < n; i++) {
      bool is_cens = (c < numc && indc(c) == i);
      if (is_cens) {
        Type Y = Type(n - i);
        for (int r = 0; r < P; r++) cum(r) += q(r, c) / (Y * Y);
        for (int r = 0; r < P; r++) psi(i, r) = q(r, c) / Y;
        c++;
      }
      for (int r = 0; r < P; r++) psi(i, r) -= cum(r);
    }

    REPORT(q);
    REPORT(psi);
  }

  return Type(0) * (beta.sum() + alpha.sum());
}
// ---- original template ends ----
} // namespace wnpmle_corr_cmprsk_BC

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
