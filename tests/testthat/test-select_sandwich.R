set.seed(42)

make_data <- function(n = 400) {
  data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
}

test_that("lm: zero-padded eta matches reduced model's fitted values", {
  dat <- make_data()
  dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(nrow(dat))

  full <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = dat)
  reduced <- lm(y ~ x1 + x2, data = dat)

  res <- select_sandwich(reduced, dat, full_formula = y ~ x1 + x2 + x3 + x4 + x5)

  expect_s3_class(res, "selectSandwich")
  expect_equal(length(res$coefficients), 6)
  expect_equal(unname(res$coefficients[c("x3", "x4", "x5")]), c(0, 0, 0))
  expect_equal(unname(res$coefficients[c("(Intercept)", "x1", "x2")]),
               unname(coef(reduced)))

  eta <- as.numeric(model.matrix(full) %*% res$coefficients)
  expect_equal(eta, unname(fitted(reduced)), tolerance = 1e-8)
})

test_that("lm: covariance matrix is symmetric, positive (semi-)definite, and unselected vars get non-zero SE", {
  dat <- make_data()
  dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(nrow(dat)) * (1 + 0.5 * abs(dat$x1))

  reduced <- lm(y ~ x1 + x2, data = dat)
  res <- select_sandwich(reduced, dat, full_formula = y ~ x1 + x2 + x3 + x4 + x5)

  V <- res$vcov
  expect_equal(V, t(V), tolerance = 1e-8)
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > -1e-8))

  expect_true(all(res$se[c("x3", "x4", "x5")] > 0))
  expect_false(res$selected["x3"])
  expect_true(res$selected["x1"])
})

test_that("lm: when nothing is dropped, result matches the classic HC0 sandwich", {
  dat <- make_data()
  dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + 0.3 * dat$x3 + rnorm(nrow(dat))

  full <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = dat)
  res <- select_sandwich(full, dat)

  X <- model.matrix(full)
  e <- residuals(full)
  bread <- solve(crossprod(X))
  meat <- crossprod(X, X * e^2)
  V_manual <- bread %*% meat %*% bread

  expect_equal(unname(res$vcov), unname(V_manual), tolerance = 1e-6)
  expect_true(all(res$selected))
})

test_that("glm (logistic): zero-padded eta matches reduced model, weights/meat are sane", {
  n <- 1000
  dat <- make_data(n)
  lin <- -0.3 + 1.2 * dat$x1 - 0.8 * dat$x2
  dat$z <- rbinom(n, 1, plogis(lin))

  full_glm <- glm(z ~ x1 + x2 + x3 + x4 + x5, data = dat, family = binomial)
  reduced_glm <- glm(z ~ x1 + x2, data = dat, family = binomial)

  res <- select_sandwich(reduced_glm, dat, full_formula = z ~ x1 + x2 + x3 + x4 + x5)

  eta <- as.numeric(model.matrix(full_glm) %*% res$coefficients)
  expect_equal(eta, unname(reduced_glm$linear.predictors), tolerance = 1e-6)

  expect_true(all(res$se > 0))
  expect_equal(res$model_type, "glm")
})

test_that("glm (logistic): when nothing is dropped, matches classic HC0 logistic sandwich", {
  n <- 1000
  dat <- make_data(n)
  lin <- -0.3 + 1.2 * dat$x1 - 0.8 * dat$x2
  dat$z <- rbinom(n, 1, plogis(lin))

  full_glm <- glm(z ~ x1 + x2 + x3 + x4 + x5, data = dat, family = binomial)
  res <- select_sandwich(full_glm, dat)

  X <- model.matrix(full_glm)
  mu <- fitted(full_glm)
  w <- mu * (1 - mu)
  s <- dat$z - mu

  A <- crossprod(X, X * w)
  B <- crossprod(X, X * s^2)
  Ainv <- solve(A)
  V_manual <- Ainv %*% B %*% Ainv

  expect_equal(unname(res$vcov), unname(V_manual), tolerance = 1e-6)
})

test_that("HC1 correction inflates variance relative to HC0", {
  dat <- make_data()
  dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(nrow(dat))
  reduced <- lm(y ~ x1 + x2, data = dat)

  res0 <- select_sandwich(reduced, dat, full_formula = y ~ x1 + x2 + x3 + x4 + x5, hc = "HC0")
  res1 <- select_sandwich(reduced, dat, full_formula = y ~ x1 + x2 + x3 + x4 + x5, hc = "HC1")

  expect_true(all(res1$se >= res0$se))
})

test_that("mismatched full_formula scope throws an informative error", {
  dat <- make_data()
  dat$y <- 1 + 2 * dat$x1 + rnorm(nrow(dat))
  reduced <- lm(y ~ x1 + x2, data = dat)  # x2 selected

  expect_error(
    select_sandwich(reduced, dat, full_formula = y ~ x1 + x3 + x4 + x5),  # x2 missing from scope
    "not found as"
  )
})

test_that("summary() and confint() run and have expected shape", {
  dat <- make_data()
  dat$y <- 1 + 2 * dat$x1 - 1.5 * dat$x2 + rnorm(nrow(dat))
  reduced <- lm(y ~ x1 + x2, data = dat)
  res <- select_sandwich(reduced, dat, full_formula = y ~ x1 + x2 + x3 + x4 + x5)

  s <- summary(res)
  expect_s3_class(s, "summary.selectSandwich")
  expect_equal(nrow(s$table), 6)

  ci <- confint(res)
  expect_equal(dim(ci), c(6, 2))
  expect_true(all(ci[, 1] <= res$coefficients + 1e-8))
  expect_true(all(ci[, 2] >= res$coefficients - 1e-8))
})
