#' Prepare reusable null-model support for TOP score scans
#'
#' Fits the nuisance model once and stores the matrices needed to scan many
#' variants. Reuse the returned object in [top_score_scan()] for successive
#' variant batches.
#'
#' @param y Numeric outcome matrix as described in [top_fit()].
#' @param baseline,additive,pairwise,saturated Numeric nuisance covariates for
#'   each second-stage structure.
#' @param missing_value Numeric missing-tumor code for cases, or `NULL`.
#' @param subtype_cutoff Subtypes with counts at or below this value are
#'   omitted.
#'
#' @return A reusable object of class `top_score_support`.
#' @export
top_score_support <- function(y, baseline = NULL, additive = NULL,
                              pairwise = NULL, saturated = NULL,
                              missing_value = NULL, subtype_cutoff = 10) {
  y <- .top_validate_outcome(y, missing_value)
  n <- nrow(y)
  baseline <- .top_covariate_matrix(baseline, "baseline", n)
  additive <- .top_covariate_matrix(additive, "additive", n)
  pairwise <- .top_covariate_matrix(pairwise, "pairwise", n)
  saturated <- .top_covariate_matrix(saturated, "saturated", n)
  .top_validate_cutoff(subtype_cutoff)

  support <- ScoreTestSupport(
    y = y,
    baselineonly = baseline,
    additive = additive,
    pairwise.interaction = pairwise,
    saturated = saturated,
    missingTumorIndicator = missing_value,
    cutoff = subtype_cutoff
  )
  required <- c(
    "inv_info_vec", "YminusP", "W_obs", "WXZ_vec", "zc", "idx.drop",
    "z.design.baselineonly", "z.design.additive",
    "z.design.pairwise.interaction", "z.design.saturated", "z.standard"
  )
  if (!is.list(support) || !all(required %in% names(support))) {
    stop("TOP could not construct a complete score-support object", call. = FALSE)
  }

  structure(
    support,
    class = c("top_score_support", "list"),
    observations = n,
    missing_value = missing_value,
    subtype_cutoff = subtype_cutoff
  )
}

#' @export
print.top_score_support <- function(x, ...) {
  cat("TOP score-scan support\n")
  cat("Observations:", attr(x, "observations"),
      " Retained subtypes:", nrow(x$z.standard), "\n")
  if (length(x$idx.drop)) {
    cat("Dropped observations:", length(x$idx.drop), "\n")
  }
  invisible(x)
}

#' Scan variants with a TOP fixed-effect score test
#'
#' Tests each variant independently using a null model that is fitted once.
#' For large scans, create `support` with [top_score_support()] and reuse it
#' across variant batches.
#'
#' @param y Numeric outcome matrix as described in [top_fit()].
#' @param variants Numeric vector, matrix, or data frame with one variant per
#'   column and one subject per row.
#' @param structure Second-stage structure used for the tested variants.
#' @param support Optional object returned by [top_score_support()].
#' @param baseline,additive,pairwise,saturated Nuisance covariates used only
#'   when `support` is not supplied.
#' @param missing_value Numeric missing-tumor code for cases, or `NULL`.
#' @param subtype_cutoff Subtypes with counts at or below this value are
#'   omitted when fitting a new support object.
#'
#' @return A `top_score_scan` object containing a tidy result table, score
#'   matrix, information array, and reusable support object.
#' @export
top_score_scan <- function(y, variants,
                           structure = c("additive", "baseline", "pairwise", "saturated"),
                           support = NULL, baseline = NULL, additive = NULL,
                           pairwise = NULL, saturated = NULL,
                           missing_value = NULL, subtype_cutoff = 10) {
  structure <- match.arg(structure)
  y <- .top_validate_outcome(y, missing_value)
  variants <- .top_variant_matrix(variants, nrow(y))

  if (is.null(support)) {
    support <- top_score_support(
      y = y,
      baseline = baseline,
      additive = additive,
      pairwise = pairwise,
      saturated = saturated,
      missing_value = missing_value,
      subtype_cutoff = subtype_cutoff
    )
  } else {
    if (!inherits(support, "top_score_support")) {
      stop("`support` must be created by `top_score_support()`", call. = FALSE)
    }
    if (!identical(attr(support, "observations"), nrow(y))) {
      stop("`support` was fitted with a different number of observations", call. = FALSE)
    }
  }

  legacy_structure <- switch(
    structure,
    baseline = "baselineonly",
    additive = "additive",
    pairwise = "pairwise.interaction",
    saturated = "saturated"
  )
  raw <- ScoreTest(
    y = y,
    x = variants,
    second.stage.structure = legacy_structure,
    score.test.support = support,
    missingTumorIndicator = missing_value
  )

  score <- as.matrix(raw$score.result)
  number_variants <- ncol(variants)
  parameters <- ncol(score)
  if (nrow(score) != number_variants ||
      !identical(dim(raw$efficient.info.result), c(number_variants * parameters, parameters))) {
    stop("TOP returned score statistics with unexpected dimensions", call. = FALSE)
  }

  design <- switch(
    structure,
    baseline = support$z.design.baselineonly,
    additive = support$z.design.additive,
    pairwise = support$z.design.pairwise.interaction,
    saturated = support$z.design.saturated
  )
  effect_names <- colnames(design)
  if (is.null(effect_names)) effect_names <- paste0("effect", seq_len(parameters))
  dimnames(score) <- list(colnames(variants), effect_names)

  information <- array(
    NA_real_,
    dim = c(parameters, parameters, number_variants),
    dimnames = list(effect_names, effect_names, colnames(variants))
  )
  statistic <- p_value <- rep(NA_real_, number_variants)
  status <- rep("ok", number_variants)
  for (index in seq_len(number_variants)) {
    rows <- ((index - 1L) * parameters + 1L):(index * parameters)
    information[, , index] <- raw$efficient.info.result[rows, , drop = FALSE]
    solved <- tryCatch(
      solve(information[, , index], score[index, ]),
      error = function(error) error
    )
    if (inherits(solved, "error")) {
      status[index] <- paste0("singular information: ", conditionMessage(solved))
    } else {
      statistic[index] <- as.numeric(crossprod(score[index, ], solved))
      p_value[index] <- stats::pchisq(statistic[index], df = parameters, lower.tail = FALSE)
    }
  }

  results <- data.frame(
    variant = colnames(variants),
    statistic = statistic,
    df = parameters,
    p.value = p_value,
    status = status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  structure(
    list(
      results = results,
      score = score,
      information = information,
      support = support,
      structure = structure
    ),
    class = "top_score_scan"
  )
}

#' @export
print.top_score_scan <- function(x, ...) {
  cat("TOP fixed-effect score scan\n")
  cat("Second-stage structure:", x$structure,
      " Variants:", nrow(x$results), "\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.top_score_scan <- function(x, ...) x$results

.top_variant_matrix <- function(variants, n) {
  if (is.data.frame(variants) && !all(vapply(variants, is.numeric, logical(1)))) {
    stop("`variants` must contain only numeric columns", call. = FALSE)
  }
  if (!is.numeric(variants)) {
    stop("`variants` must be numeric", call. = FALSE)
  }
  variants <- as.matrix(variants)
  if (nrow(variants) != n || ncol(variants) == 0L) {
    stop("`variants` must have ", n, " rows and at least one column", call. = FALSE)
  }
  if (any(!is.finite(variants))) {
    stop("`variants` cannot contain missing or non-finite values", call. = FALSE)
  }
  if (is.null(colnames(variants))) {
    colnames(variants) <- paste0("variant", seq_len(ncol(variants)))
  }
  colnames(variants) <- make.unique(colnames(variants))
  variants
}
