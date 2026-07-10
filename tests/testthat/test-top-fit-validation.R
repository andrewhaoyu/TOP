test_that("top_fit validates outcome coding before native code", {
  valid <- cbind(
    disease = c(0, 0, 1, 1),
    marker = c(NA, NA, 0, 1)
  )
  expect_error(top_fit(valid + c(0, 0, 2, 0)), "0 and 1")

  control_marker <- valid
  control_marker[1, 2] <- 0
  expect_error(top_fit(control_marker), "controls must be NA")

  case_missing <- valid
  case_missing[3, 2] <- NA
  expect_error(top_fit(case_missing), "case tumor markers cannot be NA")
})

test_that("top_fit validates covariate dimensions and cutoff", {
  y <- cbind(disease = c(0, 0, 1, 1), marker = c(NA, NA, 0, 1))
  expect_error(top_fit(y, additive = 1:3), "must have 4 rows")
  expect_error(top_fit(y, additive = c(1, 2, NA, 4)), "finite")
  expect_error(top_fit(y, subtype_cutoff = 1.5), "non-negative integer")
})
