# wnpmle 0.1.3

* The censoring-corrected sandwich variance (`se = "sandwich_adj"`, the
  default) now follows the formula in Bellach and Kosorok (2026, Appendix B.3)
  exactly. Two coding issues in the correction term were fixed: the correction
  is now subtracted from the subject-wise scores (which are on the negative
  log-likelihood scale), and two internal matrices are now filled in the right
  order. In the settings checked, the effect on the standard errors was in the
  fourth decimal place or smaller.
* Subject-wise scores and the censoring correction are now computed with TMB
  (new templates in `inst/tmb/`), which makes `se = "sandwich"` and
  `se = "sandwich_adj"` much faster for larger data sets.
* Compiled TMB templates are now kept between R sessions (in
  `tools::R_user_dir("wnpmle", "cache")`), so they are compiled only once per
  installed version instead of once per session. A different folder can be set
  with `options(wnpmle.cache_dir = ...)` or the environment variable
  `WNPMLE_CACHE_DIR`, e.g. to share one compiled copy between cluster jobs.
* The meat of the sandwich is computed with `crossprod()`, which uses much less
  memory.
* New example data `readmission_prep()`: hospital readmissions of 403
  colorectal cancer patients with death as terminal event (Gonzalez et al.,
  2005; copy of the `readmission` data from frailtypack). The vignette now uses
  these data as its main example, including model selection with
  `plot_loglik()`. The bladder quick start now compares the Ghosh-Lin fit
  with the estimates reported by Ghosh and Lin (2002).
* `predict()` now handles factor covariates (same coding as in the fit), and
  uses the correct limit `log(1 + x)` for the Box-Cox model at `rho = 0`.
* `summary()` now returns an object that is printed once (before, printing the
  result of `summary()` showed the coefficients twice).
* `plot_loglik()`: wider default PDF (`width`, `height` arguments), titles no
  longer overlap, and the points rho = 1 and r = 1 are labelled.
* Documentation: `"sandwich_adj"` is listed as the default for `se`; references now
  cite the arXiv preprint of Bellach and Kosorok (2026) with its DOI, and DOIs
  were added to all references where available.

# wnpmle 0.1.2

* First CRAN release.
