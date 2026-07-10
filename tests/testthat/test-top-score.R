test_that("score support is reusable and score scans return tidy results", {
  data("data", package = "TOP")
  y <- data[, 1:5]
  variants <- data[, 6, drop = FALSE]
  adjustment <- data[, 7, drop = FALSE]

  support <- top_score_support(
    y,
    additive = adjustment,
    missing_value = 888
  )
  scan <- top_score_scan(
    y,
    variants,
    support = support,
    missing_value = 888
  )

  expect_s3_class(support, "top_score_support")
  expect_s3_class(scan, "top_score_scan")
  expect_equal(scan$results$variant, "SNP")
  expect_equal(scan$results$df, 5)
  expect_equal(scan$results$p.value, 0.05302981, tolerance = 1e-6)
  expect_equal(as.data.frame(scan), scan$results)
  expect_equal(dim(scan$information), c(5, 5, 1))
})

test_that("complete-case score support uses the same stable fields", {
  data("data", package = "TOP")
  missing <- GenerateMissingPosition(data[, 1:5], 888)
  complete <- data[-missing, , drop = FALSE]

  support <- top_score_support(
    complete[, 1:5],
    additive = complete[, 7, drop = FALSE]
  )
  scan <- top_score_scan(
    complete[, 1:5],
    complete[, 6, drop = FALSE],
    support = support
  )

  expect_named(support)
  expect_null(support$idx.drop)
  expect_equal(scan$results$status, "ok")
  expect_true(is.finite(scan$results$p.value))
})

test_that("top_followup returns variant-specific Wald results", {
  data("data", package = "TOP")
  followup <- top_followup(
    y = data[, 1:5],
    variants = data[, 6, drop = FALSE],
    adjustment = data[, 7, drop = FALSE],
    missing_value = 888,
    keep_models = FALSE
  )

  expect_s3_class(followup, "top_followup")
  expect_equal(followup$results$variant, "SNP")
  expect_equal(followup$results$global_association_p, 0.0543, tolerance = 1e-4)
  expect_equal(followup$results$global_heterogeneity_p, 0.126, tolerance = 1e-3)
  expect_true(nrow(followup$effects) > 0)
  expect_null(followup$models)
  expect_equal(as.data.frame(followup), followup$results)
})

test_that("score and follow-up wrappers validate variant inputs", {
  y <- cbind(disease = c(0, 0, 1, 1), marker = c(NA, NA, 0, 1))
  expect_error(top_score_scan(y, variants = 1:3), "must have 4 rows")
  expect_error(top_followup(y, variants = c(1, 2, NA, 4)), "cannot contain")
})
