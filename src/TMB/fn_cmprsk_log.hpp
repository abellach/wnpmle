// fn_cmprsk_log.hpp
//
// Weighted NPMLE objective for the LOGARITHMIC competing-risks model.
// Same data layout, parameters and structure as fn_cmprsk_BC_tmb.cpp;
// only the transformation G differs:
//
//   G(x; r)  = log(1 + r x) / r      (r -> 0: G(x) = x, Fine-Gray)
//   G'(x; r) = 1 / (1 + r x)         (r = 1: proportional odds)
//
// Term-by-term (as fn() in sim_LOG_combined.R):
//   f1 = sum log(lambda)
//   f2 = sum cov1 %*% beta
//   f3 = sum log G'(exp(cov1 beta) * Lambda)          = -sum log(1 + r x)
//   f4 = sum G(exp(cova beta) * Lam.all)               over ALL n subjects
//   f5 = sum w * exp(cov2 beta) * lambda * G'(exp(cov2 beta) * Lambda)
//   nll = -f1 - f2 - f3 + f4 + f5
//
// Estimates: nlminb on (beta, alpha = log lambda).
// Inverse Fisher: solve(obj$he()) is w.r.t. (beta, alpha); convert to
// (beta, lambda) with J = diag(1, ..., 1, lambda):  J %*% solve(He) %*% J.
//

// G(x; r) = log(1 + r x) / r, with the r -> 0 limit = x
//
// Model "fn_cmprsk_log" of the wnpmle library (competing risks). Called from
// src/wnpmle.cpp when the data list passed to MakeADFun() contains
// model = "fn_cmprsk_log". The code between the two marker lines is the stand-alone
// template verified on 16 Oct 2026 (test_step12/3/4.R); only its first line
// was changed (objective_function -> model).
#ifndef FN_CMPRSK_LOG_HPP
#define FN_CMPRSK_LOG_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

namespace wnpmle_fn_cmprsk_log {
// ---- original template begins ----
template<class Type>
Type Glog(Type x, Type r) {
  if (CppAD::abs(r) < Type(1e-8)) {
    return x;
  } else {
    return log(Type(1) + r * x) / r;
  }
}

// G'(x; r) = 1 / (1 + r x)
template<class Type>
Type Glog1(Type x, Type r) {
  return Type(1) / (Type(1) + r * x);
}

template<class Type>
Type model(objective_function<Type>* obj)
{
  // ---- Data (identical to fn_cmprsk_BC_tmb.cpp, with r instead of rho) ----
  DATA_MATRIX(cova);      // n x numcov, all subjects ordered by time
  DATA_IVECTOR(idx_all);  // n, number of type-1 events at/before each subject
  DATA_MATRIX(cov1);      // num1 x numcov, type-1 subjects
  DATA_MATRIX(cov2);      // num2 x numcov, type-2 subjects
  DATA_MATRIX(w);         // num2 x num1, IPCW weights w.new
  DATA_SCALAR(r);         // log-model parameter (fixed)

  int n      = cova.rows();
  int numcov = cova.cols();
  int num1   = cov1.rows();
  int num2   = cov2.rows();

  // ---- Parameters ----
  PARAMETER_VECTOR(beta);   // numcov
  PARAMETER_VECTOR(alpha);  // num1, log(lambda)

  vector<Type> lambda(num1);
  for (int k = 0; k < num1; k++) lambda(k) = exp(alpha(k));

  vector<Type> Lambda(num1);
  Type csum = Type(0);
  for (int k = 0; k < num1; k++) { csum += lambda(k); Lambda(k) = csum; }

  // ---- f1, f2, f3: type-1 subjects ----
  Type f1 = Type(0), f2 = Type(0), f3 = Type(0);
  for (int k = 0; k < num1; k++) {
    f1 += log(lambda(k));
    Type eta1k = Type(0);
    for (int l = 0; l < numcov; l++) eta1k += cov1(k, l) * beta(l);
    f2 += eta1k;
    Type e1k = exp(eta1k);
    f3 += log(Glog1(e1k * Lambda(k), r));
  }

  // ---- f4: all n subjects ----
  Type f4 = Type(0);
  for (int i = 0; i < n; i++) {
    Type etai = Type(0);
    for (int l = 0; l < numcov; l++) etai += cova(i, l) * beta(l);
    Type ei = exp(etai);
    Type Lam_i = (idx_all(i) >= 1) ? Lambda(idx_all(i) - 1) : Type(0);
    f4 += Glog(ei * Lam_i, r);
  }

  // ---- f5: type-2 subjects x type-1 event times ----
  Type f5 = Type(0);
  for (int j = 0; j < num2; j++) {
    Type eta2j = Type(0);
    for (int l = 0; l < numcov; l++) eta2j += cov2(j, l) * beta(l);
    Type e2j = exp(eta2j);
    for (int k = 0; k < num1; k++) {
      Type wjk = w(j, k);
      if (wjk != Type(0)) {
        f5 += wjk * e2j * lambda(k) * Glog1(e2j * Lambda(k), r);
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
} // namespace wnpmle_fn_cmprsk_log

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
