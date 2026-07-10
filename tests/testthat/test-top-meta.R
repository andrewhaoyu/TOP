test_that("top_meta combines any number of studies", {
  studies <- list(
    list(score = c(1, 2), information = diag(c(2, 3))),
    list(score = c(0.5, -0.5), information = diag(c(1, 2))),
    list(score = c(-0.25, 1), information = diag(c(4, 1)))
  )
  result <- top_meta(studies)

  expect_s3_class(result, "top_meta")
  expect_equal(result$score, c(1.25, 2.5))
  expect_equal(result$information, diag(c(7, 6)))
  expect_equal(result$studies, 3)
  expect_equal(result$statistic, 1.25^2 / 7 + 2.5^2 / 6)
  expect_equal(result$p.value, pchisq(result$statistic, 2, lower.tail = FALSE))
})

test_that("two-study compatibility wrapper preserves its result", {
  old <- ScoreMetaAnalysis(c(1, 2), diag(2), c(3, 4), diag(c(2, 3)))
  expect_equal(old$score.meta, c(4, 6))
  expect_equal(old$infor.meta, diag(c(3, 4)))
})

test_that("top_meta rejects incompatible and unsafe inputs", {
  expect_error(top_meta(list()), "non-empty")
  expect_error(
    top_meta(list(list(score = 1:2, information = diag(3)))),
    "matching"
  )
  expect_error(
    top_meta(list(list(score = 1:2, information = matrix(c(1, 2, 0, 1), 2)))),
    "symmetric"
  )
  expect_error(
    top_meta(list(list(score = 1:2, information = matrix(0, 2, 2)))),
    "singular"
  )
})
