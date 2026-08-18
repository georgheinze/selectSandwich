#' Zero-padded robust sandwich covariance after variable selection
#'
#' @description
#' Takes a fitted \code{lm} or \code{glm} object that resulted from a
#' variable-selection procedure (e.g. \code{stats::step()}, manual
#' backward/forward selection, etc.) together with the data set used, and
#' returns:
#' \itemize{
#'   \item a coefficient vector over the *full* candidate set of predictors,
#'         where variables that were dropped during selection are set to
#'         exactly 0 ("zero-corrected" coefficients), and
#'   \item a robust (Huber/White sandwich-type) variance-covariance matrix
#'         for that full-length vector, in which the dropped variables get a
#'         genuine (non-zero) variance estimate rather than a deterministic
#'         zero.
#' }
#'
#' @details
#' \strong{Rationale.} Ordinary robust/sandwich standard errors are only
#' defined for the parameters that were actually estimated. If a variable
#' was dropped during selection, its "coefficient" is conventionally either
#' omitted or reported as 0 with 0 standard error -- which overstates how
#' sure we are that the true effect is zero (the variable may simply not
#' have been selected in this particular draw / this particular data set).
#'
#' \strong{Method.} Let \eqn{X} be the \eqn{n \times p} design matrix built
#' from \emph{all} candidate predictors (\code{full_formula}), and let
#' \eqn{\hat\beta_{full}} be the length-\eqn{p} vector obtained by placing
#' the fitted coefficients from the reduced (selected) model into the
#' positions of the selected variables, and 0 elsewhere. Write
#' \eqn{\hat\eta = X \hat\beta_{full}} (this is numerically identical to the
#' reduced model's linear predictor, since the zero entries do not
#' contribute). For the model family (Gaussian for \code{lm}, or the
#' \code{glm} family) define, row by row,
#' \deqn{\hat\mu_i = g^{-1}(\hat\eta_i), \quad
#'       \dot g^{-1}_i = \partial \mu / \partial \eta \big|_{\hat\eta_i}, \quad
#'       V(\hat\mu_i) \text{ the variance function.}}
#' The "bread" (expected/Fisher information built on the full design) and
#' "meat" (outer product of empirical score contributions, again built on
#' the full design) are
#' \deqn{A = X^\top W X, \qquad W = \mathrm{diag}\!\left(\dot g^{-1 2}_i /
#'       V(\hat\mu_i)\right)}
#' \deqn{B = X^\top D X, \qquad D = \mathrm{diag}\!\left((y_i-\hat\mu_i)^2
#'       \dot g^{-1 2}_i / V(\hat\mu_i)^2\right)}
#' and the zero-corrected robust covariance matrix is the usual sandwich
#' \deqn{\widehat{\mathrm{Var}}(\hat\beta_{full}) = A^{-1} B A^{-1}.}
#' For \code{lm} (identity link, constant variance function) this reduces
#' exactly to the standard Huber/White heteroskedasticity-robust covariance
#' \eqn{(X^\top X)^{-1} X^\top \mathrm{diag}(e_i^2) X (X^\top X)^{-1}}, just
#' evaluated on the full column set instead of only the selected columns.
#' For a variable that was *not* selected, its row/column of \eqn{A} and
#' \eqn{B} are non-zero whenever that variable is correlated with the
#' residuals or with the selected predictors, which is exactly why it picks
#' up a non-zero variance: the sandwich formula does not require
#' \eqn{\hat\beta_{full}} to solve the full-model score equations, only that
#' it is *some* candidate parameter vector we want the (model-robust)
#' uncertainty of.
#'
#' \strong{What you need to supply.} \code{select_sandwich()} needs to know
#' the full candidate design, i.e. the model you would have fit \emph{before}
#' selection. By default it uses every other column of \code{data} as a
#' main effect (\code{y ~ .}), which is only appropriate if that is
#' really the candidate set \code{step()} (or whatever selector you used)
#' searched over. If your selection scope included transformations,
#' interactions, or a subset of columns, pass that scope explicitly via
#' \code{full_formula}.
#'
#' @param fit A fitted model object of class \code{"lm"} or \code{"glm"}
#'   (e.g. the output of \code{stats::step()} applied to either).
#' @param data The data frame used to fit \code{fit} (or a data frame with
#'   the same columns / factor levels). Used to build the full design
#'   matrix. Rows with missing values in any variable used by
#'   \code{full_formula} are dropped (with a warning if that changes the
#'   sample size relative to \code{fit}).
#' @param full_formula Optional. A formula giving the full candidate model
#'   (the selection "scope"), e.g. \code{y ~ x1 + x2 + x3 + x4 + x5}. If
#'   \code{NULL} (default), it is built automatically as the response of
#'   \code{fit} regressed on all remaining columns of \code{data}
#'   (\code{y ~ .}).
#' @param hc Type of finite-sample correction applied to the sandwich
#'   matrix. \code{"HC0"} (default) applies none. \code{"HC1"} multiplies
#'   by \code{n / (n - k)}, where \code{k} is the number of *selected*
#'   (non-zero) coefficients -- i.e. the effective degrees of freedom used
#'   by the fitted (reduced) model.
#' @param conf_level Confidence level used by \code{confint()} on the
#'   returned object. Default 0.95.
#'
#' @return An object of class \code{"selectSandwich"}, a list with elements:
#' \describe{
#'   \item{coefficients}{Named numeric vector of length \code{p} (all
#'     candidate variables), zero for unselected variables.}
#'   \item{vcov}{The \code{p x p} zero-corrected robust covariance matrix,
#'     with matching row/column names.}
#'   \item{se}{\code{sqrt(diag(vcov))}, the zero-corrected robust standard
#'     errors (including for unselected variables).}
#'   \item{selected}{Named logical vector, \code{TRUE} for variables kept in
#'     the reduced (fitted) model.}
#'   \item{model_type}{\code{"lm"} or \code{"glm"}.}
#'   \item{family}{The (possibly synthetic, for \code{lm}) family object
#'     used.}
#'   \item{reduced_formula, full_formula}{The two formulas involved.}
#'   \item{n, p_full, p_selected}{Sample size and dimensions.}
#'   \item{hc}{The correction type used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 300
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
#'                    x4 = rnorm(n), x5 = rnorm(n))
#' dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(n)
#'
#' full  <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = dat)
#' step_fit <- step(full, direction = "backward", trace = 0)
#'
#' res <- select_sandwich(step_fit, dat)
#' print(res)
#' summary(res)
#' confint(res)
#'
#' # logistic regression example
#' dat$z <- rbinom(n, 1, plogis(-0.3 + 1.2 * dat$x1 - 0.8 * dat$x2))
#' full_glm <- glm(z ~ x1 + x2 + x3 + x4 + x5, data = dat, family = binomial)
#' step_glm <- step(full_glm, direction = "backward", trace = 0)
#' res_glm <- select_sandwich(step_glm, dat)
#' summary(res_glm)
#' }
#'
#' @export
select_sandwich <- function(fit, data, full_formula = NULL,
                             hc = c("HC0", "HC1"), conf_level = 0.95) {

  hc <- match.arg(hc)
  cl <- match.call()

  if (!inherits(fit, "lm")) {
    stop("`fit` must be an object of class 'lm' or 'glm' ",
         "(glm objects also inherit from 'lm').")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  model_type <- if (inherits(fit, "glm")) "glm" else "lm"
  family_obj <- if (model_type == "glm") fit$family else stats::gaussian()

  reduced_formula <- stats::formula(fit)
  resp_name <- .response_name(reduced_formula)

  if (is.null(full_formula)) {
    predictors <- setdiff(names(data), resp_name)
    if (length(predictors) == 0) {
      stop("Could not build a default `full_formula`: no predictor ",
           "columns found in `data` other than the response '",
           resp_name, "'. Please supply `full_formula` explicitly.")
    }
    full_formula <- stats::as.formula(
      paste(resp_name, "~", paste(predictors, collapse = " + "))
    )
  } else {
    full_formula <- stats::as.formula(full_formula)
    if (.response_name(full_formula) != resp_name) {
      stop("The response in `full_formula` ('", .response_name(full_formula),
           "') does not match the response of `fit` ('", resp_name, "').")
    }
  }

  full_terms <- stats::terms(full_formula, data = data)
  full_mf <- stats::model.frame(full_terms, data = data, na.action = stats::na.omit)

  n_fit <- length(stats::fitted(fit))
  if (nrow(full_mf) != n_fit) {
    warning("Number of complete-case rows implied by `full_formula` (",
            nrow(full_mf), ") differs from the number of observations used ",
            "to fit `fit` (", n_fit, "). Results may not be directly ",
            "comparable; make sure `data` matches what `fit` was estimated ",
            "on (e.g. no extra missingness in variables that were not ",
            "selected).")
  }

  fit_weights <- if (!is.null(fit$prior.weights)) fit$prior.weights else fit$weights
  has_case_weights <- (!is.null(stats::model.weights(full_mf))) ||
    (!is.null(fit_weights) && !all(fit_weights == 1))
  if (has_case_weights) {
    warning("Prior/case weights are not supported by select_sandwich(); ",
            "they will be ignored. Results assume unweighted estimation.")
  }

  # Reuse the reduced model's factor contrasts when available, so that
  # factor columns in the full design matrix are coded identically to how
  # `fit` coded them (otherwise coefficient names may fail to match below).
  X_full <- if (!is.null(fit$contrasts)) {
    stats::model.matrix(full_terms, data = full_mf, contrasts.arg = fit$contrasts)
  } else {
    stats::model.matrix(full_terms, data = full_mf)
  }
  y_full <- stats::model.response(full_mf)
  y_full <- .as_numeric_response(y_full)

  offset_full <- stats::model.offset(full_mf)
  if (is.null(offset_full)) offset_full <- rep(0, nrow(X_full))

  # ---- Build the zero-padded coefficient vector over the full column set --
  beta_reduced <- stats::coef(fit)
  full_names <- colnames(X_full)

  missing_names <- setdiff(names(beta_reduced), full_names)
  if (length(missing_names) > 0) {
    stop("The following coefficient name(s) from `fit` were not found as ",
         "columns of the full design matrix built from `full_formula`: ",
         paste(sQuote(missing_names), collapse = ", "), ". This usually ",
         "means `full_formula` does not actually contain the selection ",
         "scope used by `fit` (e.g. missing interaction/factor terms, or ",
         "different factor contrasts). Please supply a `full_formula` that ",
         "matches the scope that was searched over.")
  }

  beta_full <- setNames(rep(0, length(full_names)), full_names)
  beta_full[names(beta_reduced)] <- beta_reduced
  selected <- setNames(full_names %in% names(beta_reduced), full_names)

  # ---------------------------- Sandwich pieces -----------------------------
  eta <- as.numeric(X_full %*% beta_full) + offset_full

  if (model_type == "lm") {
    mu <- eta
    mu_eta <- rep(1, length(eta))
    Vmu <- rep(1, length(eta))
  } else {
    mu <- family_obj$linkinv(eta)
    mu_eta <- family_obj$mu.eta(eta)
    Vmu <- family_obj$variance(mu)
  }

  eps <- .Machine$double.eps^0.5
  Vmu[abs(Vmu) < eps] <- eps  # guard against division by ~0 near boundary

  w <- (mu_eta^2) / Vmu                       # IRLS / Fisher weights
  s <- (y_full - mu) * mu_eta / Vmu           # working score contributions

  A <- crossprod(X_full, X_full * w)          # bread:  X' W X
  B <- crossprod(X_full, X_full * s^2)        # meat:   X' D X

  A_inv <- tryCatch(solve(A), error = function(e) {
    stop("The full design matrix (from `full_formula`) is singular or ",
         "near-singular, so the sandwich 'bread' matrix cannot be inverted. ",
         "Check for perfectly collinear or constant predictors in the full ",
         "candidate set. Original error: ", conditionMessage(e))
  })

  V <- A_inv %*% B %*% A_inv

  if (hc == "HC1") {
    k <- sum(selected)
    n <- nrow(X_full)
    if (n > k) V <- V * (n / (n - k))
  }

  dimnames(V) <- list(full_names, full_names)
  V <- (V + t(V)) / 2  # enforce exact numerical symmetry

  out <- list(
    coefficients   = beta_full,
    vcov           = V,
    se             = sqrt(diag(V)),
    selected       = selected,
    model_type     = model_type,
    family         = family_obj,
    reduced_formula = reduced_formula,
    full_formula   = full_formula,
    n              = nrow(X_full),
    p_full         = length(full_names),
    p_selected     = sum(selected),
    hc             = hc,
    conf_level     = conf_level,
    call           = cl
  )
  class(out) <- "selectSandwich"
  out
}

# ------------------------------- internal helpers ---------------------------

#' @keywords internal
#' @noRd
.response_name <- function(formula) {
  tt <- stats::terms(formula)
  vars <- attr(tt, "variables")
  resp_idx <- attr(tt, "response")
  if (resp_idx == 0) {
    stop("The formula has no response variable; select_sandwich() requires ",
         "a fitted model with a response (y ~ ...).")
  }
  deparse(vars[[resp_idx + 1]])
}

#' Coerce a model.response() result to a plain numeric 0/1 (or continuous)
#' vector, mirroring what glm.fit()/lm.fit() effectively use internally.
#' @keywords internal
#' @noRd
.as_numeric_response <- function(y) {
  if (is.factor(y)) {
    if (nlevels(y) != 2) {
      stop("Factor responses with more than 2 levels are not supported. ",
           "Please convert the response to a numeric 0/1 vector before ",
           "calling select_sandwich().")
    }
    return(as.numeric(y == levels(y)[2]))
  }
  if (is.matrix(y) && ncol(y) == 2) {
    # cbind(success, failure) style binomial response
    total <- rowSums(y)
    total[total == 0] <- NA
    return(y[, 1] / total)
  }
  if (is.logical(y)) return(as.numeric(y))
  as.numeric(y)
}
