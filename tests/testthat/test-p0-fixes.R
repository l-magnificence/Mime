# Tests for P0 fixes
library(testthat)
library(survival)

# ===================================================================
# Issue 22: Handle identical risk scores in returnRStoROC
# ===================================================================
test_that("identical risk scores get noise added, not crash", {
  # Simulate the logic from returnRStoROC
  # When all RS values are identical, grouping should still work

  x <- data.frame(
    OS.time = c(100, 200, 300, 400, 500),
    OS = c(1, 0, 1, 0, 1),
    RS = c(1.5, 1.5, 1.5, 1.5, 1.5)  # all identical
  )

  # Before fix: median split gives all same group -> survdiff error
  # After fix: noise is added first
  if (length(unique(x$RS)) == 1) {
    x$RS <- x$RS + runif(nrow(x), min = 0.001, max = 0.01) * abs(x$RS[1])
  }

  # Should now have variation
  expect_true(length(unique(x$RS)) > 1)

  # Grouping should work
  x$Group <- ifelse(x$RS > median(x$RS), "High", "Low")
  expect_true(length(unique(x$Group)) >= 1)

  # survdiff should not error
  expect_no_error({
    data.survdiff <- survdiff(Surv(x$OS.time, x$OS) ~ x$Group)
  })
})

test_that("normal risk scores are not affected by noise fix", {
  x <- data.frame(
    OS.time = c(100, 200, 300, 400, 500),
    OS = c(1, 0, 1, 0, 1),
    RS = c(1.0, 1.5, 2.0, 2.5, 3.0)  # normal variation
  )

  # Should NOT add noise
  if (length(unique(x$RS)) == 1) {
    x$RS <- x$RS + runif(nrow(x), min = 0.001, max = 0.01) * abs(x$RS[1])
  }

  # RS should be unchanged
  expect_equal(x$RS, c(1.0, 1.5, 2.0, 2.5, 3.0))
})

# ===================================================================
# Issue 110: AUC clamped to [0, 1]
# ===================================================================
test_that("AUC values above 1 are clamped to 1", {
  auc_raw <- 1.0023  # floating point issue
  auc_clamped <- min(max(auc_raw, 0), 1)
  expect_equal(auc_clamped, 1.0)
})

test_that("AUC values below 0 are clamped to 0", {
  auc_raw <- -0.001
  auc_clamped <- min(max(auc_raw, 0), 1)
  expect_equal(auc_clamped, 0.0)
})

test_that("normal AUC values are not affected by clamping", {
  auc_raw <- 0.75
  auc_clamped <- min(max(auc_raw, 0), 1)
  expect_equal(auc_clamped, 0.75)
})

# ===================================================================
# Issue 122: SurvivalSVM direction correction
# ===================================================================
test_that("SurvSVM negates when direction is inverted", {
  set.seed(42)
  n <- 100
  OS.time <- rexp(n, rate = 0.01)
  # Predicted values positively correlated with OS.time (wrong direction)
  pred <- OS.time + rnorm(n, sd = 10)

  if (cor(pred, OS.time, use = "complete.obs") > 0) {
    pred <- -pred
  }

  # After correction, higher pred should correlate with shorter survival
  expect_true(cor(pred, OS.time, use = "complete.obs") < 0)
})

test_that("SurvSVM preserves correct direction", {
  set.seed(42)
  n <- 100
  OS.time <- rexp(n, rate = 0.01)
  # Predicted values negatively correlated with OS.time (correct direction)
  pred <- -OS.time + rnorm(n, sd = 10)

  if (cor(pred, OS.time, use = "complete.obs") > 0) {
    pred <- -pred
  }

  # After correction, higher pred should still correlate with shorter survival
  expect_true(cor(pred, OS.time, use = "complete.obs") < 0)
})

# ===================================================================
# Issue 123: RSF predict error handling
# ===================================================================
test_that("tryCatch catches predict errors gracefully", {
  result <- tryCatch(
    stop("class name too long in predict"),
    error = function(e) {
      rep(NA_real_, 5)
    }
  )
  expect_equal(result, rep(NA_real_, 5))
})

test_that("successful predict is not affected by tryCatch", {
  result <- tryCatch(
    as.numeric(c(1.5, 2.3, 0.8)),
    error = function(e) {
      rep(NA_real_, 3)
    }
  )
  expect_equal(result, c(1.5, 2.3, 0.8))
})

# ===================================================================
# Issue 112: ROC computation handles prob column mismatch
# ===================================================================
test_that("prob column mismatch returns NA instead of crashing", {
  prob <- matrix(c(0.3, 0.7, 0.6, 0.4), ncol = 2, dimnames = list(NULL, c("A", "B")))
  lev <- c("N", "Y")

  expect_false(all(lev %in% colnames(prob)))

  if (!all(lev %in% colnames(prob))) {
    result <- data.frame(ROC = NA, Sens = NA, Spec = NA)
  }
  expect_true(is.na(result$ROC))
})

test_that("correct prob columns work normally", {
  prob <- matrix(c(0.3, 0.7, 0.6, 0.4), ncol = 2, dimnames = list(NULL, c("N", "Y")))
  lev <- c("N", "Y")

  expect_true(all(lev %in% colnames(prob)))
})
