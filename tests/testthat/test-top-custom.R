test_that("custom subtype designs return a top_model", {
  data("data", package = "TOP")
  y <- data[, 1:5]
  z_standard <- GenerateZstandard(y, missingTumorIndicator = 888)
  design <- matrix(0, nrow(z_standard), 5)
  colnames(design) <- c("low", "positive", "high", "her2", "negative")

  hr_positive <- z_standard[, 1] == 1 | z_standard[, 2] == 1
  her2_positive <- z_standard[, 3] == 1
  high_grade <- z_standard[, 4] == 2
  design[hr_positive & !her2_positive & !high_grade, 1] <- 1
  design[hr_positive & her2_positive, 2] <- 1
  design[hr_positive & !her2_positive & high_grade, 3] <- 1
  design[!hr_positive & her2_positive, 4] <- 1
  design[!hr_positive & !her2_positive, 5] <- 1

  fit <- top_fit_custom(
    y,
    covariate = data[, 6, drop = FALSE],
    design = design,
    adjustment = data[, 7, drop = FALSE],
    missing_value = 888
  )

  expect_s3_class(fit, "top_model")
  expect_equal(attr(fit, "custom_design"), design)
  expect_equal(nrow(fit$global.test), 2)
  expect_equal(unique(fit$second.stage.test$Covariate), c("SNP", "PC1"))
  expect_equal(unique(fit$first.stage.test$Covariate), c("SNP", "PC1"))
})

test_that("custom subtype design dimensions are validated", {
  y <- cbind(disease = c(0, 0, 1, 1), marker = c(NA, NA, 0, 1))
  expect_error(
    top_fit_custom(y, 1:4, matrix(1, 3, 1), subtype_cutoff = 0),
    "must have"
  )
})
