#define TMB_LIB_INIT R_init_fn_BC_tmb
#include <TMB.hpp>

template<class Type>
inline Type log1p_ad(const Type &x) {
  return log(Type(1) + x);
}

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

  const int num1 = cov1.rows();
  const int num2 = cov2.rows();
  const int n02  = cov02.rows();

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

  vector<Type> lambda = exp(alpha);
  vector<Type> Lambda(num1);
  Type csum = Type(0);
  for(int j=0; j<num1; ++j){ csum += lambda(j); Lambda(j) = csum; }

  vector<Type> eta1  = cov1  * beta;
  vector<Type> eta2  = cov2  * beta;
  vector<Type> eta02 = cov02 * beta;
  vector<Type> e1  = exp(eta1);
  vector<Type> e2  = exp(eta2);
  vector<Type> e02 = exp(eta02);

  Type f1 = alpha.sum();
  Type f2 = eta1.sum();

  Type f3 = Type(0);
  for(int j=0; j<num1; ++j) f3 += log1p_ad(e1(j) * Lambda(j));
  f3 *= (rho - Type(1));

  Type f4 = Type(0);
  for(int u=0; u<n02; ++u){
    int k = idx02(u);
    Type Lam4u = (k >= 0) ? Lambda(k) : Type(0);
    Type base  = Type(1) + e02(u) * Lam4u;
    if (CppAD::abs(rho) < Type(1e-10)) { f4 += log(base); }
    else { f4 += (pow(base, rho) - Type(1)) / rho; }
  }

  Type f5 = Type(0);
  for(int i=0; i<num2; ++i){
    Type vi  = e2(i);
    Type den = kmc2(i);
    int start = idx2(i) + 1;
    if(start < 0) start = 0;
    if(start >= num1) continue;
    for(int j=start; j<num1; ++j){
      Type w = kmc1(j) / den;
      Type oneplus = Type(1) + vi * Lambda(j);
      f5 += w * vi * lambda(j) * pow(oneplus, rho - Type(1));
    }
  }

  Type nll = -f1 - f2 - f3 + f4 + f5;
  REPORT(f1); REPORT(f2); REPORT(f3); REPORT(f4); REPORT(f5); REPORT(nll);
  ADREPORT(beta); ADREPORT(alpha); ADREPORT(lambda); ADREPORT(Lambda);
  return nll;
}
