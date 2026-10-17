# TMB models in wnpmle

Each model is one file in this folder. All of them are compiled into the
package's library when the package is installed; `src/wnpmle.cpp` only
connects the model names to the files. In R, a model is used as

```r
TMB::MakeADFun(data = list(model = "<name>", ...), parameters, DLL = "wnpmle")
```

and the name is chosen by `.wnpmle_tmb_model()` in `R/helpers.R` (recurrent
events) or `.wnpmle_tmb_model_cmprsk()` in `R/cmprsk_fit.R` (competing risks).

| Name | File | Data | Model | Computes | Used by |
|---|---|---|---|---|---|
| `fn_BC` | `fn_BC.hpp` | recurrent events | Box-Cox (`rho`) | negative log-likelihood; estimates and inverse Fisher | `wnpmle_fit()` |
| `fn_log` | `fn_log.hpp` | recurrent events | logarithmic (`r`) | negative log-likelihood; estimates and inverse Fisher | `wnpmle_fit()` |
| `score_BC` | `score_BC.hpp` | recurrent events | Box-Cox | subject-wise scores (plain sandwich) | `.compute_score_tmb()` |
| `score_log` | `score_log.hpp` | recurrent events | logarithmic | subject-wise scores (plain sandwich) | `.compute_score_tmb()` |
| `corr_BC` | `corr_BC.hpp` | recurrent events | Box-Cox | censoring correction psi (`sandwich_adj`) | `.censoring_correction_tmb()` |
| `corr_log` | `corr_log.hpp` | recurrent events | logarithmic | censoring correction psi (`sandwich_adj`) | `.censoring_correction_tmb()` |
| `fn_cmprsk_BC` | `fn_cmprsk_BC.hpp` | competing risks | Box-Cox (`rho`) | negative log-likelihood; estimates and inverse Fisher | `.wnpmle_fit_cmprsk()` |
| `fn_cmprsk_log` | `fn_cmprsk_log.hpp` | competing risks | logarithmic (`r`) | negative log-likelihood; estimates and inverse Fisher | `.wnpmle_fit_cmprsk()` |
| `score_cmprsk_BC` | `score_cmprsk_BC.hpp` | competing risks | Box-Cox | subject-wise scores (plain sandwich) | `.wnpmle_fit_cmprsk()` |
| `score_cmprsk_log` | `score_cmprsk_log.hpp` | competing risks | logarithmic | subject-wise scores (plain sandwich) | `.wnpmle_fit_cmprsk()` |
| `corr_cmprsk_BC` | `corr_cmprsk_BC.hpp` | competing risks | Box-Cox | censoring correction psi (`sandwich_adj`) | `.wnpmle_fit_cmprsk()` |
| `corr_cmprsk_log` | `corr_cmprsk_log.hpp` | competing risks | logarithmic | censoring correction psi (`sandwich_adj`) | `.wnpmle_fit_cmprsk()` |

## How a model file is built

Each file is a stand-alone TMB template with a short wrapper around it:

```cpp
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
namespace wnpmle_<name> {
// ---- original template begins ----
template<class Type>
Type model(objective_function<Type>* obj)   // was: objective_function<Type>::operator()()
{ ... unchanged ... }
// ---- original template ends ----
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
```

The namespace keeps helper functions of different models apart, so two files
may use the same helper name.

## Adding a model

1. Put it in a new file `<name>.hpp` here, with the wrapper above.
2. In `src/wnpmle.cpp`, add `#include "TMB/<name>.hpp"` and
   `if (model == "<name>") return wnpmle_<name>::model(this);`.
3. Add a row to the table above, and a test in `tests/testthat/`.

The existing files do not need to change.
