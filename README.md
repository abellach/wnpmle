# wnpmle

**Weighted NPMLE for Recurrent Events with a Competing Terminal Event**

An R package implementing the weighted nonparametric maximum likelihood
estimator (wNPMLE) for the marginal mean of recurrent events when a competing
terminal event (e.g. death) is present.

## Background

In many clinical studies, subjects experience repeated events (e.g. COPD
exacerbations, bladder tumour recurrences) until a terminal event permanently
ends the process. Treating death as independent censoring leads to
overestimation of the expected number of recurrences. This package implements
the semiparametric estimator of Bellach and Kosorok (2026), which correctly
accounts for the terminal event through inverse probability of censoring
weights and a weighted likelihood.

Two transformation models are provided:

- **Box-Cox**: G(x) = ((1 + x)^rho - 1) / rho — reduces to the **Ghosh-Lin model** when rho = 1
- **Logarithmic**: G(x) = log(1 + r*x) / r — reduces to the **proportional odds model** when r = 1

When no terminal events are present, the package automatically reduces to the
plain NPMLE, equivalent to the **Zeng-Lin model** (Zeng and Lin, 2006).

Both models are estimated via automatic differentiation (TMB). Standard errors
are available via the sandwich estimator with optional censoring correction.

## Installation

```r
# install.packages("remotes")
remotes::install_github("abellach/wnpmle")
```

**Requirements**: The package requires [TMB](https://github.com/kaskr/adcomp)
and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Windows only).

## Quick start

```r
library(wnpmle)

# Prepare bladder cancer data
bdata       <- bladder_prep()
bdata_clean <- bdata[, c("id", "time", "status", "treat", "num", "size")]

# --- Ghosh-Lin model (Box-Cox, rho = 1) ---
fit_bc <- wnpmle_fit(
  Surv(time, status) ~ treat + num + size,
  data  = bdata_clean,
  id    = "id",
  model = "boxcox",
  rho   = 1,
  se    = "sandwich"
)
summary(fit_bc)
plot(fit_bc)
baseline(fit_bc)
AIC(fit_bc)
BIC(fit_bc)

# --- Proportional odds model (Log, r = 1) ---
fit_log <- wnpmle_fit(
  Surv(time, status) ~ treat + num + size,
  data  = bdata_clean,
  id    = "id",
  model = "log",
  rho   = 1,
  se    = "sandwich"
)
summary(fit_log)
plot(fit_log)
baseline(fit_log)
AIC(fit_log)
BIC(fit_log)

# --- Profile log-likelihood for transformation parameter ---
result <- plot_loglik(
  Surv(time, status) ~ treat + num + size,
  data     = bdata_clean,
  id       = "id",
  tau      = 59,
  rho_grid = seq(0.01, 1.2, by = 0.01),
  r_grid   = seq(0.01, 1.2, by = 0.01)
)

# Table of log-likelihood values
head(result)

# Optimal transformation parameters
result[result$model == "boxcox", ][which.max(result$loglik[result$model == "boxcox"]), ]
result[result$model == "log",    ][which.max(result$loglik[result$model == "log"   ]), ]
```

## Status codes

| Status | Meaning                     |
|--------|-----------------------------|
| 0      | Censored                    |
| 1      | Recurrent event             |
| 2      | Terminal event (e.g. death) |

## Data format

Users can bring their own data. The required format is one row per event per
subject with columns:

- `id` — subject identifier
- `time` — event time
- `status` — 0 (censored), 1 (recurrent event), 2 (terminal event)
- covariates of interest

## Citation

If you use this package in your research, please cite:

Bellach, A. and Kosorok, M.R. (2026). Weighted NPMLE for the marginal mean of
recurrent events with a competing terminal event. 

https://arxiv.org/abs/2605.25934

### BibTeX

```bibtex
@article{bellach2026wnpmle,
  title   = {Weighted {NPMLE} for the marginal mean of recurrent events
             with a competing terminal event},
  author  = {Bellach, Anna and Kosorok, Michael R.},
  year   = {2026},
  url    = {https://arxiv.org/abs/2605.25934}
}
```

## References

Bellach A., Kosorok M.R., Rüschendorf L. and Fine J.P. (2019). Weighted NPMLE for the
subdistribution of a competing risk. *JASA*, 114(525), 259-270.

Ghosh, D. and Lin, D.Y. (2002). Marginal regression models for recurrent and
terminal events. *Statistica Sinica*, 12, 663-688.

Zeng, D. and Lin, D.Y. (2006). Efficient estimation of semiparametric
transformation models for counting processes. *Biometrika*, 93(3), 627-640.

## License

GPL (>= 3)

