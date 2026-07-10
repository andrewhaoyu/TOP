test_that("the documented breast example fits end to end", {
  data("data", package = "TOP")

  fit <- top_fit(
    y = data[, 1:5],
    additive = data[, 6:7],
    missing_value = 888,
    subtype_cutoff = 10
  )

  expect_s3_class(fit, "top_model")
  expect_named(fit)
  expect_length(coef(fit), 31)
  expect_equal(dim(vcov(fit)), c(31, 31))
  expect_equal(as.numeric(fit$loglikelihood), -16726.47, tolerance = 0.05)
  expect_equal(as.numeric(fit$AIC), 33514.95, tolerance = 0.05)
  expect_s3_class(summary(fit), "summary.top_model")
})
