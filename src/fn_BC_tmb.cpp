#include <TMB.hpp>

/*
 Box-Cox transformation model (TMB)
 Parameters: beta (p), alpha (num1) with lambda = exp(alpha), Lambda = cumsum(lambda)
 Data: cov1, cov2, cov02, idx02, idx2, kmc1, kmc2, rho
 
 Index conventions expected:
 idx02 in 0..num1   (R: findInterval(ind02, ind1))
 idx2  in -1..num1-1 (R: findInterval(ind2, ind1) - 1)
 */

template<class Type>
inline Type log1p_ad(const Type &x) {
  // AD-safe log(1 + x)
  return log(Type(1) + x);
}

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
  
  DATA_SCALAR(rho);     // fixed for profiling / grid
  
  const int num1 = cov1.rows();
  const int num2 = cov2.rows();
  const int n02  = cov02.rows();
  
  PARAMETER_VECTOR(beta);   // length p
  PARAMETER_VECTOR(alpha);  // length num1
  
  // baseline increments and cumulative baseline
  vector<Type> lambda = exp(alpha);
  vector<Type> Lambda(num1);
  Type csum = Type(0);
  for(int j=0; j<num1; ++j){
    csum += lambda(j);
    Lambda(j) = csum;
  }
  
  // linear predictors and exp()
  vector<Type> eta1  = cov1  * beta;
  vector<Type> eta2  = cov2  * beta;
  vector<Type> eta02 = cov02 * beta;
  
  vector<Type> e1  = exp(eta1);
  vector<Type> e2  = exp(eta2);
  vector<Type> e02 = exp(eta02);
  
  // f1 = sum(log(lambda)) = sum(alpha)
  Type f1 = alpha.sum();
  
  // f2 = sum(cov1 %*% beta)
  Type f2 = eta1.sum();
  
  // f3 = (rho - 1) * sum_j log(1 + exp(eta1_j) * Lambda_j)
  Type f3 = Type(0);
  for(int j=0; j<num1; ++j){
    f3 += log1p_ad(e1(j) * Lambda(j));
  }
  f3 *= (rho - Type(1));

  // f4 = sum_u ((1 + exp(eta02_u) * Lam4_u)^rho - 1)/rho
  // Lam4_u = Lambda[idx02[u]] if idx02[u] >= 0 else 0
  // idx02 in R: findInterval(ind02, ind1) - 1L  (values -1..num1-1)
  Type f4 = Type(0);
  for(int u=0; u<n02; ++u){
    int k = idx02(u);
    Type Lam4u = (k >= 0) ? Lambda(k) : Type(0);
    Type base  = Type(1) + e02(u) * Lam4u;
    
    // stable limit as rho -> 0: (base^rho - 1)/rho -> log(base)
    if (CppAD::abs(rho) < Type(1e-10)) {
      f4 += log(base);
    } else {
      f4 += (pow(base, rho) - Type(1)) / rho;
    }
  }  
  
  // f5 = sum_i sum_{j > idx2[i]} (kmc1[j]/kmc2[i]) * exp(eta2[i]) * lambda[j] * (1 + exp(eta2[i]) * Lambda[j])^(rho-1)
  Type f5 = Type(0);
  for(int i=0; i<num2; ++i){
    Type vi  = e2(i);
    Type den = kmc2(i);
    
    int start = idx2(i) + 1;
    if(start < 0) start = 0;
    if(start >= num1) continue;
    
    for(int j=start; j<num1; ++j){
      Type w = kmc1(j) / den;
      Type oneplus = Type(1) + vi * Lambda(j);   // <-- key change
      f5 += w * vi * lambda(j) * pow(oneplus, rho - Type(1));
    }
  }
  
  // nll = -f1 - f2 - f3 + f4 + f5
  Type nll = -f1 - f2 - f3 + f4 + f5;
  
  // Debug reports so you can match R term-by-term
  REPORT(f1); REPORT(f2); REPORT(f3); REPORT(f4); REPORT(f5);
  REPORT(nll);
  
  ADREPORT(beta);
  ADREPORT(alpha);
  ADREPORT(lambda);
  ADREPORT(Lambda);
  
  return nll;
}

