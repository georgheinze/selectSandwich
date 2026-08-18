\name{select_sandwich}
\alias{select_sandwich}
\title{Zero-padded robust sandwich covariance after variable selection}
\usage{
select_sandwich(fit, data, full_formula = NULL,
                 hc = c("HC0", "HC1"), conf_level = 0.95)
}
\arguments{
\item{fit}{A fitted model object of class \code{"lm"} or \code{"glm"}
  (e.g. the output of \code{stats::step()} applied to either).}

\item{data}{The data frame used to fit \code{fit}. Used to build the full
  design matrix. Rows with missing values in any variable used by
  \code{full_formula} are dropped (with a warning if that changes the
  sample size relative to \code{fit}).}

\item{full_formula}{Optional. A formula giving the full candidate model
  (the selection "scope"), e.g. \code{y ~ x1 + x2 + x3 + x4 + x5}. If
  \code{NULL} (default), it is built automatically as the response of
  \code{fit} regressed on all remaining columns of \code{data}
  (\code{y ~ .}).}

\item{hc}{Type of finite-sample correction applied to the sandwich matrix.
  \code{"HC0"} (default) applies none. \code{"HC1"} multiplies by
  \code{n / (n - k)}, where \code{k} is the number of selected (non-zero)
  coefficients.}

\item{conf_level}{Confidence level used by \code{confint()} on the returned
  object. Default 0.95.}
}
\value{
An object of class \code{"selectSandwich"}, a list with elements
\code{coefficients} (full-length, zero-padded coefficient vector),
\code{vcov} (the corresponding robust covariance matrix), \code{se},
\code{selected} (logical, which variables were kept by the selection
procedure), \code{model_type} (\code{"lm"} or \code{"glm"}), \code{family},
\code{reduced_formula}, \code{full_formula}, \code{n}, \code{p_full},
\code{p_selected}, \code{hc}, and \code{call}.
}
\description{
Given a fitted \code{lm}/\code{glm} object obtained after variable
selection (e.g. via \code{stats::step()}), returns a coefficient vector
over the \emph{full} candidate predictor set (0 for unselected variables)
together with a robust sandwich-type variance-covariance matrix computed
on the full design. Because the sandwich "bread" and "meat" are evaluated
using the full set of candidate columns, unselected variables receive a
genuine, non-zero variance estimate reflecting the fact that their
non-selection is not certain proof of a null effect.
}
\details{
See the package README / source for the full mathematical definition.
In brief: let \eqn{X} be the design matrix built from all candidate
predictors, and \eqn{\hat\beta} the length-\eqn{p} vector with the
reduced model's fitted coefficients in the selected positions and 0
elsewhere. The bread is \eqn{A = X'WX} and the meat is \eqn{B = X'DX},
with \eqn{W} and \eqn{D} built from the family's link/variance functions
evaluated at \eqn{X\hat\beta}; the returned covariance is
\eqn{A^{-1}BA^{-1}}. For \code{lm} this is exactly the usual
heteroskedasticity-robust (Huber/White) sandwich, just extended to the
unselected columns.
}
\examples{
\dontrun{
set.seed(1)
n <- 300
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                   x4 = rnorm(n), x5 = rnorm(n))
dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(n)

full <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = dat)
step_fit <- step(full, direction = "backward", trace = 0)

res <- select_sandwich(step_fit, dat)
summary(res)
confint(res)
}
}
