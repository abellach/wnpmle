# wnpmle

**Weighted NPMLE for Recurrent Events with a Competing Terminal Event**

An R package implementing the weighted nonparametric maximum likelihood
estimator (wNPMLE) for the marginal mean of recurrent events when a competing
terminal event (e.g. death) is present.

## Background

In many clinical settings, subjects experience repeated events (e.g. COPD
exacerbations, bladder tumour recurrences) until a terminal event permanently
ends the process. Treating death as independent censoring leads to
overestimation of the expected number of recurrences. This package implements
the semiparametric estimator of Bellach and Kosorok (2026), which correctly
accounts for the terminal event through inverse probability of censoring
weights and a weighted likelihood.

Two transformation models are provided:

- **Box-Cox**: G(x) = ((1 + x)^ρ − 1) / ρ
- **Logarithmic**: G(x) = log(1 + r·x) / r

Both are estimated via automatic differentiation (TMB). Standard errors are
available via the sandwich estimator with optional censoring correction.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("bellach/wnpmle")
```

**Requirements**: The package links to [TMB](https://github.com/kaskr/adcomp).
Install it first:

```r
install.packages("TMB")
```

## Quick start

```r
library(wnpmle)

# Prepare bladder cancer data
bdata <- bladder_prep()

# Fit the logarithmic transformation model
fit <- wnpmle_fit(
  Surv(time, status) ~ rx + size + number,
  data  = bdata,
  id    = "id",
  model = "log",
  rho   = 1,
  se    = "sandwich"
)

summary(fit)
plot(fit)

# Fit the Box-Cox transformation model
fit_bc <- wnpmle_fit(
  Surv(time, status) ~ rx + size + number,
  data  = bdata,
  id    = "id",
  model = "boxcox",
  rho   = 0.5,
  se    = "sandwich_corrected"
)

summary(fit_bc)
```

## Status codes

| Status | Meaning                    |
|--------|----------------------------|
| 0      | Censored                   |
| 1      | Recurrent event            |
| 2      | Terminal event (e.g. death)|

## References

Bellach, A. and Kosorok, M.R. (2026). Weighted NPMLE for the marginal mean of
recurrent events with a competing terminal event. *Journal of the American
Statistical Association*, to appear.

Bellach, A., Kosorok, M.R., Fine, J.P. (2019). Weighted NPMLE for the
subdistribution of a competing risk. *JASA*, 114(525), 259–270.

## License

GPL (≥ 3)
