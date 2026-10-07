// score_BC_tmb.cpp
//
// Subject-wise likelihood contributions for the Box-Cox recurrent-event
// model, for the SANDWICH variance only (the fit still uses fn_BC_tmb.cpp,
// which is not touched).
//
// Translation of the three f.ind blocks in bladder_BoxCox_feb2016.R
// (lines 220-270). Each subject's contribution is collected in nll_i:
//
//   censored row u (status 0):   G(e02_u * Lambda at u)                         [f.ind, lines 221-232]
//   terminal row i (status 2):   G(e2_i * Lambda at i)                          [f.ind, lines 240-252]
//                              + sum_j w_ij e2_i lambda_j (1 + e2_i Lambda_j)^(rho-1)
//   recurrent row j (status 1): -log(lambda_j) - eta1_j
//                              - (rho-1) log(1 + e1_j Lambda_j)                 [f.ind, lines 257-268]
//
// G(x) = ((1+x)^rho - 1)/rho, with G(x) = log(1+x) for rho -> 0.
//
// Usage in R: MakeADFun(..., ADreport = TRUE). Then obj$gr(par) is the
// numi x (numcov + num1) matrix of subject-wise derivatives with respect to
// (beta, alpha), alpha = log(lambda). Divide the alpha columns by lambda to
// get derivatives with respect to lambda, as in gradi.
//
// Data are the same as for fn_BC_tmb.cpp, plus the subject index of each row:
//   id1  (num1): subject of each recurrent-event row, 0-based (M1$id - 1)
//   id02 (n02):  subject of each censored/terminal row, 0-based (M02$id - 1)
//   id2  (num2): subject of each terminal row, 0-based (M2$id - 1)
//   numi:        number of subjects
//
#define TMB_LIB_INIT R_init_score_BC_tmb
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

  // ---- recurrent-event rows: -f1 - f2 - f3 ----
  for (int j = 0; j < num1; ++j) {
    nll_i(id1(j)) += -alpha(j) - eta1(j)
                     - (rho - Type(1)) * log(Type(1) + e1(j) * Lambda(j));
  }

  // ---- censored and terminal rows: G term (f4) ----
  // same indexing as fn_BC_tmb.cpp: idx02 in -1..num1-1
  for (int u = 0; u < n02; ++u) {
    int  k     = idx02(u);
    Type Lam4u = (k >= 0) ? Lambda(k) : Type(0);
    Type base  = Type(1) + e02(u) * Lam4u;
    Type G;
    if (CppAD::abs(rho) < Type(1e-10)) { G = log(base); }
    else { G = (pow(base, rho) - Type(1)) / rho; }
    nll_i(id02(u)) += G;
  }

  // ---- terminal rows: weighted term (f5) ----
  // same indexing as fn_BC_tmb.cpp: recurrent events after the terminal row
  for (int i = 0; i < num2; ++i) {
    Type vi  = e2(i);
    Type den = kmc2(i);
    int start = idx2(i) + 1;
    if (start < 0) start = 0;
    if (start >= num1) continue;
    for (int j = start; j < num1; ++j) {
      Type w       = kmc1(j) / den;
      Type oneplus = Type(1) + vi * Lambda(j);
      nll_i(id2(i)) += w * vi * lambda(j) * pow(oneplus, rho - Type(1));
    }
  }

  Type nll = nll_i.sum();   // equals nll of fn_BC_tmb.cpp (check in R)

  REPORT(nll_i);
  REPORT(nll);
  ADREPORT(nll_i);
  return nll;
}
