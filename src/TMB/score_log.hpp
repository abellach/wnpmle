// score_log_tmb.cpp
//
// Subject-wise likelihood contributions for the logarithmic recurrent-event
// model, for the SANDWICH variance only (the fit still uses fn_log_tmb.cpp,
// which is not touched).
//
// Translation of the three f.ind blocks in bladder_logarithmic_mar2016.R
// (lines 217-273). Each subject's contribution is collected in nll_i:
//
//   censored row u (status 0):   log(1 + r e02_u Lambda at u) / r                [f.ind, lines 218-229]
//   terminal row i (status 2):   log(1 + r e2_i Lambda at i) / r                 [f.ind, lines 236-253]
//                              + sum_j w_ij e2_i lambda_j / (1 + r e2_i Lambda_j)
//   recurrent row j (status 1): -log(lambda_j) - eta1_j
//                              + log(1 + r e1_j Lambda_j)                        [f.ind, lines 258-271]
//
// Usage in R: MakeADFun(..., ADreport = TRUE). Then obj$gr(par) is the
// numi x (numcov + num1) matrix of subject-wise derivatives with respect to
// (beta, alpha), alpha = log(lambda). Divide the alpha columns by lambda to
// get derivatives with respect to lambda, as in gradi.
//
// Data are the same as for fn_log_tmb.cpp, plus the subject index of each row:
//   id1  (num1): subject of each recurrent-event row, 0-based (M1$id - 1)
//   id02 (n02):  subject of each censored/terminal row, 0-based (M02$id - 1)
//   id2  (num2): subject of each terminal row, 0-based (M2$id - 1)
//   numi:        number of subjects
//
#define TMB_LIB_INIT R_init_score_log_tmb
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
  DATA_SCALAR(r);
  DATA_IVECTOR(id1);
  DATA_IVECTOR(id02);
  DATA_IVECTOR(id2);
  DATA_INTEGER(numi);

  const int num1 = cov1.rows();
  const int num2 = cov2.rows();
  const int n02  = cov02.rows();

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

  vector<Type> lambda = exp(alpha);
  vector<Type> Lambda(num1);
  Type csum = Type(0);
  for (int j = 0; j < num1; ++j) { csum += lambda(j); Lambda(j) = csum; }

  vector<Type> eta1  = cov1  * beta;
  vector<Type> eta2  = cov2  * beta;
  vector<Type> eta02 = cov02 * beta;
  vector<Type> e1  = exp(eta1);
  vector<Type> e2  = exp(eta2);
  vector<Type> e02 = exp(eta02);

  vector<Type> nll_i(numi);
  nll_i.setZero();

  // ---- recurrent-event rows: -f1 - f2 - f3, with f3 = -log(1 + r e1 Lambda) ----
  for (int j = 0; j < num1; ++j) {
    nll_i(id1(j)) += -alpha(j) - eta1(j)
                     + log(Type(1) + r * e1(j) * Lambda(j));
  }

  // ---- censored and terminal rows: G term (f4) ----
  // same indexing as fn_log_tmb.cpp: idx02 in 0..num1
  for (int u = 0; u < n02; ++u) {
    int  k     = idx02(u);
    Type Lam4u = (k >= 1) ? Lambda(k - 1) : Type(0);
    nll_i(id02(u)) += log(Type(1) + r * e02(u) * Lam4u) / r;
  }

  // ---- terminal rows: weighted term (f5) ----
  for (int i = 0; i < num2; ++i) {
    Type vi  = e2(i);
    Type den = kmc2(i);
    int start = idx2(i) + 1;
    if (start < 0) start = 0;
    if (start >= num1) continue;
    for (int j = start; j < num1; ++j) {
      Type w       = kmc1(j) / den;
      Type oneplus = Type(1) + r * vi * Lambda(j);
      nll_i(id2(i)) += w * vi * lambda(j) / oneplus;
    }
  }

  Type nll = nll_i.sum();   // equals nll of fn_log_tmb.cpp (check in R)

  REPORT(nll_i);
  REPORT(nll);
  ADREPORT(nll_i);
  return nll;
}
