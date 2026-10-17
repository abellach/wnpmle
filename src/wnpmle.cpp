// wnpmle.cpp -- the single compiled library of the wnpmle package.
//
// Each model lives in its own file in src/TMB/ (see src/TMB/README.md).
// This file only connects model names to those files: R code calls
//   TMB::MakeADFun(data = list(model = "fn_BC", ...), ..., DLL = "wnpmle")
// and the matching model below is run.
//
// To add a model: put it in src/TMB/<name>.hpp (same wrapper as the others),
// add one #include line and one "if" line below, and one row in the README.

#define TMB_LIB_INIT R_init_wnpmle
#include <TMB.hpp>

// recurrent events with a competing terminal event
#include "TMB/fn_BC.hpp"
#include "TMB/fn_log.hpp"
#include "TMB/score_BC.hpp"
#include "TMB/score_log.hpp"
#include "TMB/corr_BC.hpp"
#include "TMB/corr_log.hpp"

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_STRING(model);
  if (model == "fn_BC")     return wnpmle_fn_BC::model(this);
  if (model == "fn_log")    return wnpmle_fn_log::model(this);
  if (model == "score_BC")  return wnpmle_score_BC::model(this);
  if (model == "score_log") return wnpmle_score_log::model(this);
  if (model == "corr_BC")   return wnpmle_corr_BC::model(this);
  if (model == "corr_log")  return wnpmle_corr_log::model(this);
  error("wnpmle: unknown TMB model");
  return Type(0);
}
