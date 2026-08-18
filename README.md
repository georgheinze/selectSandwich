# selectSandwich

Robust (sandwich-type) variance-covariance estimation **after variable
selection**, for `lm` and `glm` models -- e.g. after running `stats::step()`.

## The problem

After `step()` (or any other selection procedure) drops some variables, the
fitted model object only has coefficients -- and standard errors -- for the
*selected* variables. Reporting "0, SE = 0" for the dropped variables
overstates how sure we are that their true effect is exactly zero: they were
simply not selected in *this* fit, not proven to be null.

## The idea

`select_sandwich()`:

1. Builds the design matrix for the **full candidate model** (all variables
   that were eligible for selection, selected or not).
2. Zero-pads the fitted coefficients into that full-length vector (selected
   variables keep their estimate; unselected variables get exactly 0).
3. Computes a Huber/White **sandwich** covariance matrix using the bread and
   meat evaluated on the **full** design matrix at that zero-padded vector,
   rather than on only the selected columns.

Because the bread/meat use the full covariate set, unselected variables end
up with a genuine (small, but typically non-zero) variance -- driven by how
correlated they are with the model's residuals and with the selected
predictors -- instead of a deterministic zero. For linear models this
reduces exactly to the ordinary heteroskedasticity-robust (HC0) sandwich
covariance, just extended to the dropped columns; for GLMs it uses the
family's link/variance functions in the usual IRLS-weighted sandwich form.

This is a *model-robust* variance estimate at a fixed, chosen parameter
vector (the zero-padded one) -- it does not require that vector to solve the
full model's score equations, only that we want its (heteroskedasticity- and
misspecification-robust) uncertainty.

## Installation

```r
# from a local clone / the extracted folder:
devtools::install("path/to/selectSandwich")
# or, without devtools:
install.packages("path/to/selectSandwich", repos = NULL, type = "source")
```

The package has no dependencies beyond base R (`stats`); `testthat` is only
needed to run the included unit tests.

## Usage

```r
library(selectSandwich)

set.seed(1)
n <- 300
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                   x4 = rnorm(n), x5 = rnorm(n))
dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(n)

full_fit <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = dat)
step_fit <- step(full_fit, direction = "backward", trace = 0)

res <- select_sandwich(step_fit, dat)
res            # zero-corrected coefficients
summary(res)   # Estimate / robust SE / t or z / p-value, for ALL candidates
confint(res)   # confidence intervals, including for unselected variables
```

`select_sandwich()` auto-detects `lm` vs `glm` from the class of `fit`, so
the same call works for a logistic regression selected with `step()`:

```r
dat$z <- rbinom(n, 1, plogis(-0.3 + 1.2 * dat$x1 - 0.8 * dat$x2))
full_glm <- glm(z ~ x1 + x2 + x3 + x4 + x5, data = dat, family = binomial)
step_glm <- step(full_glm, direction = "backward", trace = 0)

res_glm <- select_sandwich(step_glm, dat)
summary(res_glm)
```

### Specifying the full candidate scope

By default, `select_sandwich()` assumes the full candidate model was
`y ~ .` (every other column of `data` as a main effect). If the actual
selection scope was different (a subset of columns, interactions,
transformations, `poly()`, etc.), pass it explicitly:

```r
res <- select_sandwich(
  step_fit, dat,
  full_formula = y ~ x1 + x2 + x3 + x4 + x5 + x1:x2
)
```

If the coefficient names in `fit` can't be matched to columns of the full
design (e.g. `full_formula` doesn't actually cover the scope that was
searched, or factor contrasts differ), `select_sandwich()` stops with an
explicit error naming the mismatched terms rather than silently producing
wrong output.

### HC0 vs HC1

`hc = "HC0"` (default) applies no finite-sample correction. `hc = "HC1"`
multiplies the sandwich matrix by `n / (n - k)`, where `k` is the number of
*selected* (non-zero) coefficients -- the effective degrees of freedom used
by the fitted model.

## What this does *not* handle

- Case/prior weights (`weights = ` in `lm`/`glm`) are ignored, with a
  warning, if present.
- Factor responses with more than 2 levels, or multi-column (`cbind`)
  binomial responses, are not supported.
- This only accounts for uncertainty *given* the full candidate set you
  specify; it does not attempt to model the full stochastic selection
  procedure itself (e.g. repeated-sampling selection frequency).

## Files

```
DESCRIPTION, NAMESPACE, LICENSE
R/select_sandwich.R   -- master function + internal helpers
R/methods.R            -- print / summary / coef / vcov / confint methods
man/select_sandwich.Rd -- documentation
tests/testthat/        -- unit tests (incl. cross-checks against the
                           classic HC0 sandwich when no variables are
                           dropped)
```
