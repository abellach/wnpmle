#define TMB_LIB_INIT R_init_fn_log_tmb
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

  Type f1 = alpha.sum();
  Type f2 = eta1.sum();

  Type f3 = Type(0);
  for (int j = 0; j < num1; ++j) f3 += log(Type(1) + r * e1(j) * Lambda(j));
  f3 *= Type(-1);

  Type f4 = Type(0);
  for (int u = 0; u < n02; ++u) {
    int  k     = idx02(u);
    Type Lam4u = (k >= 1) ? Lambda(k - 1) : Type(0);
    Type arg   = r * e02(u) * Lam4u;
    f4 += log(Type(1) + arg) / r;
  }

  Type f5 = Type(0);
  for (int i = 0; i < num2; ++i) {
    Type vi  = e2(i);
    Type den = kmc2(i);
    int start = idx2(i) + 1;
    if (start < 0)    start = 0;
    if (start >= num1) continue;
    for (int j = start; j < num1; ++j) {
      Type w       = kmc1(j) / den;
      Type oneplus = Type(1) + r * vi * Lambda(j);
      f5 += w * vi * lambda(j) / oneplus;
    }
  }

  Type nll = -f1 - f2 - f3 + f4 + f5;
  REPORT(f1); REPORT(f2); REPORT(f3); REPORT(f4); REPORT(f5); REPORT(nll);
  ADREPORT(beta); ADREPORT(alpha); ADREPORT(lambda); ADREPORT(Lambda);
  return nll;
}
