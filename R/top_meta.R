#' Combine TOP score statistics across studies
#'
#' Combines independent study-level score vectors and information matrices.
#' Any number of studies may be supplied, and study order does not affect the
#' result.
#'
#' @param studies A non-empty list. Each element must be a list containing a
#'   numeric `score` vector and a square numeric `information` matrix.
#'
#' @return An object of class `top_meta` containing the combined score,
#'   information matrix, chi-squared statistic, degrees of freedom, p-value,
#'   and number of studies.
#' @export
#'
#' @examples
#' studies <- list(
#'   list(score = c(1, 2), information = diag(c(2, 3))),
#'   list(score = c(0.5, -0.5), information = diag(c(1, 2)))
#' )
#' top_meta(studies)
top_meta <- function(studies) {
  if (!is.list(studies) || length(studies) == 0L) {
    stop("`studies` must be a non-empty list", call. = FALSE)
  }

  validated <- lapply(seq_along(studies), function(index) {
    study <- studies[[index]]
    label <- paste0("studies[[", index, "]] ")
    if (!is.list(study) || !all(c("score", "information") %in% names(study))) {
      stop(label, "must contain `score` and `information`", call. = FALSE)
    }

    score <- study$score
    information <- study$information
    if (!is.numeric(score) || !is.null(dim(score)) || length(score) == 0L ||
        any(!is.finite(score))) {
      stop(label, "`score` must be a finite, non-empty numeric vector", call. = FALSE)
    }
    if (!is.numeric(information) || length(dim(information)) != 2L ||
        nrow(information) != length(score) || ncol(information) != length(score) ||
        any(!is.finite(information))) {
      stop(label, "`information` must be a finite square matrix matching `score`", call. = FALSE)
    }
    if (!isTRUE(all.equal(information, t(information), tolerance = 1e-8,
                          check.attributes = FALSE))) {
      stop(label, "`information` must be symmetric", call. = FALSE)
    }
    list(score = as.numeric(score), information = as.matrix(information))
  })

  dimensions <- vapply(validated, function(study) length(study$score), integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("all studies must use score vectors with the same length", call. = FALSE)
  }

  score <- Reduce(`+`, lapply(validated, `[[`, "score"))
  information <- Reduce(`+`, lapply(validated, `[[`, "information"))
  solution <- tryCatch(
    solve(information, score),
    error = function(error) {
      stop("combined information matrix is singular: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  statistic <- as.numeric(crossprod(score, solution))
  degrees_freedom <- length(score)

  structure(
    list(
      score = score,
      information = information,
      statistic = statistic,
      df = degrees_freedom,
      p.value = stats::pchisq(statistic, df = degrees_freedom, lower.tail = FALSE),
      studies = length(validated)
    ),
    class = "top_meta"
  )
}

#' @export
print.top_meta <- function(x, ...) {
  cat("TOP score meta-analysis\n")
  cat("Studies:", x$studies, "\n")
  cat("Chi-squared:", format(x$statistic, digits = 5),
      "on", x$df, "df; p =", format.pval(x$p.value, digits = 4), "\n")
  invisible(x)
}
