test_that("registered dosage conversion is available", {
  probabilities <- c(
    0.7, 0.2, 0.1,
    0.1, 0.3, 0.6
  )
  expect_equal(TOP:::convert(probabilities, 2), c(0.4, 1.5))
})

test_that("native routines are registered and dynamic lookup is disabled", {
  routines <- getDLLRegisteredRoutines("TOP")$.C
  expect_true(all(c("Mvpoly_complete", "EMStep", "ScoreTest") %in% names(routines)))
  expect_false(getLoadedDLLs()[["TOP"]][["dynamicLookup"]])
})
