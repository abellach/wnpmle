#include <TMB.hpp>

/*
 Logarithmic transformation model (TMB)
 Transformation: G(x) = log(1 + r*x) / r
 Survival:       S(t|Z) = (1 + r * exp(eta) * Lambda(t))^{-1/r}

 Parameters: beta (p), alpha (num1) with lambda = exp(alpha), Lambda = cumsum(lambda)
 Data: cov1, cov2, cov02, idx02, idx2, kmc1, kmc2, r

 Index conventions (same as BC code):
   idx02 in 0..num1   (R: findInterval(ind02, ind1))
   idx2  in -1..num1-1 (R: findInterval(ind2, ind1) - 1L)

 Five log-likelihood contributions:
   f1 = sum(alpha)                                              [type-1 jump sizes]
   f2 = sum(eta1)                                               [covariate term at type-1 events]
   f3 = -sum_j log(1 + r * e^{eta1_j} * Lambda_j)              [log G' at type-1 events]
   f4 = sum_u log(1 + r * e^{eta02_u} * Lambda_{idx02_u}) / r  [G(e^eta02 * Lambda) term]
   f5 = sum_i sum_{j > idx2_i} (kmc1_j / kmc2_i) *
            e^{eta2_i} * lambda_j / (1 + r * e^{eta2_i} * Lambda_j)   [type-2/cens term]

 nll = -f1 - f2 - f3 + f4 + f5
*/

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(cov1);    // num1 x p
  DATA_MATRIX(cov2);    // num2 x p
  DATA_MATRIX(cov02);   // n02  x p

  DATA_IVECTOR(idx02);  // length n02, values 0..num1
  DATA_IVECTOR(idx2);   // length num2, values -1..num1-1

  DATA_VECTOR(kmc1);    // length num1
  DATA_VECTOR(kmc2);    // length num2

  DATA_SCALAR(r);       // logarithmic transformation parameter (fixed)

  const int num1 = cov1.rows();
  const int num2 = cov2.rows();
  const int n02  = cov02.rows();

  PARAMETER_VECTOR(beta);   // length p
  PARAMETER_VECTOR(alpha);  // length num1

  // --- baseline increments and cumulative baseline ---
  vector<Type> lambda = exp(alpha);
  vector<Type> Lambda(num1);
  Type csum = Type(0);
  for (int j = 0; j < num1; ++j) {
    csum   += lambda(j);
    Lambda(j) = csum;
  }

  // --- linear predictors ---
  vector<Type> eta1  = cov1  * beta;
  vector<Type> eta2  = cov2  * beta;
  vector<Type> eta02 = cov02 * beta;

  vector<Type> e1  = exp(eta1);
  vector<Type> e2  = exp(eta2);
  vector<Type> e02 = exp(eta02);

  // -------------------------------------------------------------------
  // f1 = sum(log lambda_j) = sum(alpha_j)
  // -------------------------------------------------------------------
  Type f1 = alpha.sum();

  // -------------------------------------------------------------------
  // f2 = sum_j eta1_j   (covariate contribution at type-1 event times)
  // -------------------------------------------------------------------
  Type f2 = eta1.sum();

  // -------------------------------------------------------------------
  // f3 = -sum_j log(1 + r * e^{eta1_j} * Lambda_j)
  //
  // G'(x) = 1 / (1 + r*x)
  //   => log G'(e1_j * Lambda_j) = -log(1 + r * e1_j * Lambda_j)
  // Lambda_j = cumsum including jump j (right-continuous),
  // consistent with the reference fn and the R-side sandwich.
  // -------------------------------------------------------------------
  Type f3 = Type(0);
  for (int j = 0; j < num1; ++j) {
    f3 += log(Type(1) + r * e1(j) * Lambda(j));
  }
  f3 *= Type(-1);

  // -------------------------------------------------------------------
  // f4 = (1/r) * sum_u log(1 + r * e^{eta02_u} * Lambda_{idx02_u})
  //
  // This is -log S(X_u | Z_u) = G(e^{eta02_u} * Lambda_{idx02_u})
  //   = log(1 + r * e^{eta02_u} * Lambda_{idx02_u}) / r
  //
  // idx02(u) = k:  if k == 0, no type-1 event before u => Lambda = 0
  //                if k >= 1, Lambda = Lambda[k-1]  (0-based: Lambda[k-1])
  //
  // Convention (same as BC):
  //   idx02 comes from R's findInterval(ind02, ind1), giving values 0..num1.
  //   k == 0  => no type-1 event precedes u  => Lam4u = 0
  //   k >= 1  => Lam4u = Lambda[k-1]  (Lambda is 0-indexed in C++)
  // -------------------------------------------------------------------
  Type f4 = Type(0);
  for (int u = 0; u < n02; ++u) {
    int  k     = idx02(u);
    Type Lam4u = (k >= 1) ? Lambda(k - 1) : Type(0);
    Type arg   = r * e02(u) * Lam4u;
    f4 += log(Type(1) + arg) / r;
  }

  // -------------------------------------------------------------------
  // f5 = sum_i sum_{j > idx2_i} (kmc1_j / kmc2_i) *
  //              e^{eta2_i} * lambda_j / (1 + r * e^{eta2_i} * Lambda_j)
  //
  // This is the IPCW-weighted integral of dG(e^{eta2_i} * Lambda(t)) for
  // type-2 / censored subjects.  G'(x) = 1/(1+rx), so the integrand is
  //   e^{eta2_i} * lambda_j * G'(e^{eta2_i} * Lambda_j)
  //   = e^{eta2_i} * lambda_j / (1 + r * e^{eta2_i} * Lambda_j)
  //
  // idx2(i) in -1..num1-1; start from idx2(i)+1 (0-based j).
  // -------------------------------------------------------------------
  Type f5 = Type(0);
  for (int i = 0; i < num2; ++i) {
    Type vi  = e2(i);
    Type den = kmc2(i);

    int start = idx2(i) + 1;
    if (start < 0)    start = 0;
    if (start >= num1) continue;

    for (int j = start; j < num1; ++j) {
      Type w        = kmc1(j) / den;
      Type oneplus  = Type(1) + r * vi * Lambda(j);
      f5 += w * vi * lambda(j) / oneplus;
    }
  }

  // -------------------------------------------------------------------
  // nll = -f1 - f2 - f3 + f4 + f5
  //   (f3 is already negative, so -f3 adds its absolute value)
  // -------------------------------------------------------------------
  Type nll = -f1 - f2 - f3 + f4 + f5;

  // --- Debug reports (match R term-by-term) ---
  REPORT(f1); REPORT(f2); REPORT(f3); REPORT(f4); REPORT(f5);
  REPORT(nll);

  ADREPORT(beta);
  ADREPORT(alpha);
  ADREPORT(lambda);
  ADREPORT(Lambda);

  return nll;
}
