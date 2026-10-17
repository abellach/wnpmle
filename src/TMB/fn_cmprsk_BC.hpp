// fn_cmprsk_BC.hpp
//
// TMB translation of fn() from sim_BC_combined2.R (lines 162-181).
// Weighted NPMLE objective for the Box-Cox competing-risks model.
// rho = 1 -> Fine-Gray; rho -> 0 -> proportional odds.
//
// Term-by-term correspondence with the R fn():
//   f1 <- sum(log(lambda))                                      -> f1 below
//   f2 <- sum(cov1 %*% beta)                                    -> f2 below
//   f3 <- sum(log((1+exp(cov1%*%beta)*Lambda)^{rho-1}))         -> f3 below
//   f4 <- sum(((1+exp(cova%*%beta)*Lam.all)^rho-1)/rho)         -> f4 below (loops over ALL n subjects)
//   f5 <- sum(w.new*V5*Lamb*(1+V5*Lam)^{rho-1})                 -> f5 below
//   f  <- -f1-f2-f3+f4+f5
//
// AD gives the gradient and Hessian automatically, so this replaces both
// optim(fn, gr=NULL, method="BFGS") AND the separate hessian(fback, ...)
// call in the R script -- sdreport() below returns the Fisher-inverse SEs
// directly (this is "var.bread" in the R script).
//

// G(x; rho) = ((1+x)^rho - 1)/rho, with the rho -> 0 limit = log(1+x)
//
// Model "fn_cmprsk_BC" of the wnpmle library (competing risks). Called from
// src/wnpmle.cpp when the data list passed to MakeADFun() contains
// model = "fn_cmprsk_BC". The code between the two marker lines is the stand-alone
// template verified on 16 Oct 2026 (test_step12/3/4.R); only its first line
// was changed (objective_function -> model).
#ifndef FN_CMPRSK_BC_HPP
#define FN_CMPRSK_BC_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

namespace wnpmle_fn_cmprsk_BC {
// ---- original template begins ----
template<class Type>
Type Gbc(Type x, Type rho) {
  if (CppAD::abs(rho) < Type(1e-8)) {
    return log(Type(1) + x);
  } else {
    return (pow(Type(1) + x, rho) - Type(1)) / rho;
  }
}

template<class Type>
Type model(objective_function<Type>* obj)
{
  // ---- Data ----
  // cova: n x numcov, covariates for ALL subjects, ordered by time (== M$z1,z2 in R)
  // idx_all: length n, 0-indexed count of type-1 events with time <= subject's own
  //          time (0 = "no type-1 event has occurred yet"); this is what lets us
  //          look up Lam.all[i] = Lambda[idx_all[i]-1] without needing the M4 matrix.
  // cov1: num1 x numcov, covariates of the type-1 (event of interest) subjects
  // cov2: num2 x numcov, covariates of the type-2 (competing event) subjects
  // w:    num2 x num1, the IPCW weight matrix w.new from the R script
  // rho:  Box-Cox transformation parameter (fixed, not estimated)
  DATA_MATRIX(cova);
  DATA_IVECTOR(idx_all);
  DATA_MATRIX(cov1);
  DATA_MATRIX(cov2);
  DATA_MATRIX(w);
  DATA_SCALAR(rho);

  int n      = cova.rows();
  int numcov = cova.cols();
  int num1   = cov1.rows();
  int num2   = cov2.rows();

  // ---- Parameters ----
  PARAMETER_VECTOR(beta);   // length numcov
  PARAMETER_VECTOR(alpha);  // length num1, log(lambda) -- keeps lambda > 0 automatically

  vector<Type> lambda(num1);
  for (int k = 0; k < num1; k++) lambda(k) = exp(alpha(k));

  // Lambda(k) = cumulative jump size up through the k-th type-1 event time
  // (this is Ma3 %*% lambda in the R script, but Ma3 is just a lower-triangular
  // ones matrix, so a running cumulative sum is exactly equivalent and O(num1)
  // instead of O(num1^2)).
  vector<Type> Lambda(num1);
  Type csum = Type(0);
  for (int k = 0; k < num1; k++) { csum += lambda(k); Lambda(k) = csum; }

  // ---- f1, f2, f3: loop over type-1 subjects ----
  Type f1 = Type(0), f2 = Type(0), f3 = Type(0);
  for (int k = 0; k < num1; k++) {
    f1 += log(lambda(k));

    Type eta1k = Type(0);
    for (int l = 0; l < numcov; l++) eta1k += cov1(k, l) * beta(l);
    f2 += eta1k;

    Type e1k = exp(eta1k);
    f3 += (rho - Type(1)) * log(Type(1) + e1k * Lambda(k));
  }

  // ---- f4: loop over ALL n subjects (matches sum over cova/Lam.all in R) ----
  Type f4 = Type(0);
  for (int i = 0; i < n; i++) {
    Type etai = Type(0);
    for (int l = 0; l < numcov; l++) etai += cova(i, l) * beta(l);
    Type ei = exp(etai);
    Type Lam_i = (idx_all(i) >= 1) ? Lambda(idx_all(i) - 1) : Type(0);
    f4 += Gbc(ei * Lam_i, rho);
  }

  // ---- f5: double loop over (type-2 subject j) x (type-1 event k) ----
  Type f5 = Type(0);
  for (int j = 0; j < num2; j++) {
    Type eta2j = Type(0);
    for (int l = 0; l < numcov; l++) eta2j += cov2(j, l) * beta(l);
    Type e2j = exp(eta2j);
    for (int k = 0; k < num1; k++) {
      Type wjk = w(j, k);
      if (wjk != Type(0)) {
        f5 += wjk * e2j * lambda(k) * pow(Type(1) + e2j * Lambda(k), rho - Type(1));
      }
    }
  }

  Type nll = -f1 - f2 - f3 + f4 + f5;

  ADREPORT(beta);
  ADREPORT(lambda);
  ADREPORT(Lambda);

  return nll;
}
// ---- original template ends ----
} // namespace wnpmle_fn_cmprsk_BC

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
