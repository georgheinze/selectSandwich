#' @export
coef.selectSandwich <- function(object, ...) object$coefficients

#' @export
vcov.selectSandwich <- function(object, ...) object$vcov

#' @export
print.selectSandwich <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Zero-corrected coefficients after variable selection (",
      x$model_type, ", ", x$p_selected, "/", x$p_full,
      " variables selected)\n", sep = "")
  cat("Robust covariance type:", x$hc, "\n\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' Summary table for a selectSandwich object
#'
#' Produces the usual Estimate / Std. Error / test statistic / p-value
#' table, computed for *every* candidate variable (including those that
#' were not selected, which get the zero-corrected robust standard error
#' rather than an omitted row or a deterministic 0).
#'
#' @param object A \code{"selectSandwich"} object.
#' @param ... Unused, for S3 consistency.
#' @export
summary.selectSandwich <- function(object, ...) {
  est <- object$coefficients
  se <- object$se
  use_t <- object$model_type == "lm"

  stat <- est / se
  if (use_t) {
    df_resid <- object$n - object$p_selected
    pval <- 2 * stats::pt(-abs(stat), df = max(df_resid, 1))
    stat_name <- "t value"
  } else {
    pval <- 2 * stats::pnorm(-abs(stat))
    stat_name <- "z value"
  }

  tab <- cbind(
    Estimate   = est,
    `Std. Error` = se,
    stat,
    `Pr(>|.|)` = pval
  )
  colnames(tab)[3] <- stat_name
  rownames(tab) <- names(est)

  out <- list(
    table      = tab,
    selected   = object$selected,
    model_type = object$model_type,
    hc         = object$hc,
    n          = object$n,
    p_full     = object$p_full,
    p_selected = object$p_selected,
    reduced_formula = object$reduced_formula,
    full_formula    = object$full_formula,
    call       = object$call
  )
  class(out) <- "summary.selectSandwich"
  out
}

#' @export
print.summary.selectSandwich <- function(x, digits = max(3L, getOption("digits") - 3L),
                                          signif.stars = getOption("show.signif.stars"), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nModel type:", x$model_type, " | Robust covariance:", x$hc, "\n")
  cat("Selected", x$p_selected, "of", x$p_full, "candidate variables; n =", x$n, "\n\n")

  tab <- x$table
  rn <- rownames(tab)
  rownames(tab) <- ifelse(x$selected[rn], rn, paste0(rn, " (unselected)"))

  stats::printCoefmat(tab, digits = digits, signif.stars = signif.stars,
                       na.print = "NA", has.Pvalue = TRUE)

  cat("\nNote: rows marked '(unselected)' were dropped during variable\n",
      "selection and have coefficient 0 by construction; their standard\n",
      "error reflects sampling/selection uncertainty about that zero,\n",
      "not the precision of an estimated effect.\n", sep = "")
  invisible(x)
}

#' Confidence intervals for a selectSandwich object
#'
#' @param object A \code{"selectSandwich"} object.
#' @param parm Which parameters to compute intervals for (names or indices);
#'   defaults to all.
#' @param level Confidence level; defaults to \code{object$conf_level}.
#' @param ... Unused.
#' @export
confint.selectSandwich <- function(object, parm, level = object$conf_level, ...) {
  est <- object$coefficients
  se <- object$se
  if (missing(parm)) parm <- names(est)

  use_t <- object$model_type == "lm"
  alpha <- 1 - level
  if (use_t) {
    df_resid <- max(object$n - object$p_selected, 1)
    crit <- stats::qt(1 - alpha / 2, df = df_resid)
  } else {
    crit <- stats::qnorm(1 - alpha / 2)
  }

  lo <- est[parm] - crit * se[parm]
  hi <- est[parm] + crit * se[parm]
  out <- cbind(lo, hi)
  pct <- paste0(format(100 * c(alpha / 2, 1 - alpha / 2), trim = TRUE, digits = 3), " %")
  colnames(out) <- pct
  rownames(out) <- parm
  out
}
