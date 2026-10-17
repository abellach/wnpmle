// score_cmprsk_BC.hpp
//
// Subject-wise contributions to the negative log-likelihood of the Box-Cox
// competing-risks model (same likelihood as fn_cmprsk_BC_tmb.cpp).
// Used with MakeADFun(..., ADreport = TRUE): obj$gr() then returns the
// n x (numcov + num1) matrix of subject-wise scores w.r.t. (beta, alpha);
// dividing the alpha columns by lambda gives the scores w.r.t. lambda.
//
// Contribution of subject i (time-ordered position i, eta_i = Z_i beta,
// x_i = exp(eta_i) * Lambda(X_i), Lambda(X_i) = sum of jumps up to X_i):
//   all subjects:     + G(x_i)
//   type 1 (k-th):    - log lambda_k - eta_i - log G'(x_i)
//   type 2 (j-th):    + sum_k w(j,k) exp(eta_i) lambda_k G'(exp(eta_i) Lambda_k)
// with G(x) = ((1+x)^rho - 1)/rho, log G'(x) = (rho - 1) log(1 + x).
// The contributions sum to the objective in fn_cmprsk_BC_tmb.cpp.
//
//
// Model "score_cmprsk_BC" of the wnpmle library (competing risks). Called from
// src/wnpmle.cpp when the data list passed to MakeADFun() contains
// model = "score_cmprsk_BC". The code between the two marker lines is the stand-alone
// template verified on 16 Oct 2026 (test_step12/3/4.R); only its first line
// was changed (objective_function -> model).
#ifndef SCORE_CMPRSK_BC_HPP
#define SCORE_CMPRSK_BC_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

namespace wnpmle_score_cmprsk_BC {
// ---- original template begins ----
template<class Type>
Type Gbc(Type x, Type rho) {
  if (CppAD::abs(rho) < Type(1e-8)) return log(Type(1) + x);
  return (pow(Type(1) + x, rho) - Type(1)) / rho;
}

template<class Type>
Type model(objective_function<Type>* obj)
{
  DATA_MATRIX(cova);      // n x numcov, all subjects ordered by time
  DATA_IVECTOR(idx_all);  // n, number of type-1 events at/before each subject
  DATA_IVECTOR(status);   // n, 0 / 1 / 2
  DATA_IVECTOR(j2);       // n, row in w for type-2 subjects (0-based), -1 otherwise
  DATA_MATRIX(w);         // num2 x num1, IPCW weights
  DATA_SCALAR(rho);

  int n = cova.rows(), numcov = cova.cols(), num1 = w.cols();

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

  vector<Type> lambda(num1), Lambda(num1);
  Type csum = Type(0);
  for (int k = 0; k < num1; k++) { lambda(k) = exp(alpha(k)); csum += lambda(k); Lambda(k) = csum; }

  vector<Type> contrib(n);
  for (int i = 0; i < n; i++) {
    Type eta = Type(0);
    for (int l = 0; l < numcov; l++) eta += cova(i, l) * beta(l);
    Type e  = exp(eta);
    int  m  = idx_all(i);
    Type Li = (m >= 1) ? Lambda(m - 1) : Type(0);
    Type x  = e * Li;

    Type c = Gbc(x, rho);
    if (status(i) == 1) {
      c += -log(lambda(m - 1)) - eta - (rho - Type(1)) * log(Type(1) + x);
    }
    if (status(i) == 2) {
      int j = j2(i);
      for (int k = 0; k < num1; k++) {
        if (w(j, k) != Type(0))
          c += w(j, k) * e * lambda(k) * pow(Type(1) + e * Lambda(k), rho - Type(1));
      }
    }
    contrib(i) = c;
  }

  ADREPORT(contrib);
  return contrib.sum();
}
// ---- original template ends ----
} // namespace wnpmle_score_cmprsk_BC

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
