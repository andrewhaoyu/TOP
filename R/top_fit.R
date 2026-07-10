.top_covariate_matrix <- function(value, name, n) {
  if (is.null(value)) {
    return(NULL)
  }
  if (is.data.frame(value) && !all(vapply(value, is.numeric, logical(1)))) {
    stop("`", name, "` must contain only numeric covariates", call. = FALSE)
  }
  if (!is.numeric(value)) {
    stop("`", name, "` must be numeric", call. = FALSE)
  }
  value <- as.matrix(value)
  if (nrow(value) != n) {
    stop("`", name, "` must have ", n, " rows", call. = FALSE)
  }
  if (ncol(value) == 0L || any(!is.finite(value))) {
    stop("`", name, "` must contain finite numeric covariates", call. = FALSE)
  }
  if (is.null(colnames(value))) {
    colnames(value) <- paste0(name, seq_len(ncol(value)))
  }
  value
}

.top_validate_outcome <- function(y, missing_value) {
  if (is.data.frame(y) && !all(vapply(y, is.numeric, logical(1)))) {
    stop("`y` must contain only numeric columns", call. = FALSE)
  }
  if (!is.numeric(y) || length(dim(y)) != 2L || ncol(y) < 2L || nrow(y) < 2L) {
    stop("`y` must be a numeric matrix with a disease column and at least one tumor marker", call. = FALSE)
  }
  y <- as.matrix(y)
  disease <- y[, 1L]
  if (anyNA(disease) || !all(disease %in% c(0, 1))) {
    stop("the first column of `y` must contain only case-control values 0 and 1", call. = FALSE)
  }
  if (length(unique(disease)) != 2L) {
    stop("`y` must contain at least one case and one control", call. = FALSE)
  }
  tumor <- y[, -1L, drop = FALSE]
  if (any(!is.na(tumor[disease == 0, , drop = FALSE]))) {
    stop("tumor markers for controls must be NA", call. = FALSE)
  }
  if (is.null(missing_value) && anyNA(tumor[disease == 1, , drop = FALSE])) {
    stop("case tumor markers cannot be NA unless a missing-value analysis is specified", call. = FALSE)
  }
  if (!is.null(missing_value) &&
      (!is.numeric(missing_value) || length(missing_value) != 1L || !is.finite(missing_value))) {
    stop("`missing_value` must be NULL or one finite numeric value", call. = FALSE)
  }
  y
}

#' Fit a TOP model with validated inputs and named results
#'
#' This is the recommended interface for new analyses. The original
#' [TwoStageModel()] function remains available for backward compatibility.
#'
#' @param y Numeric matrix or data frame whose first column is case-control
#'   status (0/1) and remaining columns are tumor characteristics. Tumor
#'   characteristics must be `NA` for controls.
#' @param baseline Covariates with a common baseline effect across subtypes.
#' @param additive Covariates using the additive second-stage model.
#' @param pairwise Covariates using the pairwise-interaction model.
#' @param saturated Covariates using the saturated model.
#' @param missing_value Numeric code used for missing tumor characteristics in
#'   cases, or `NULL` when case tumor data are complete.
#' @param start Optional starting vector for second-stage parameters.
#' @param subtype_cutoff Subtypes with counts at or below this value are
#'   omitted, matching the historical TOP behavior.
#'
#' @return A named list of class `top_model`.
#' @export
top_fit <- function(y, baseline = NULL, additive = NULL, pairwise = NULL,
                    saturated = NULL, missing_value = NULL, start = NULL,
                    subtype_cutoff = 10) {
  y <- .top_validate_outcome(y, missing_value)
  n <- nrow(y)
  baseline <- .top_covariate_matrix(baseline, "baseline", n)
  additive <- .top_covariate_matrix(additive, "additive", n)
  pairwise <- .top_covariate_matrix(pairwise, "pairwise", n)
  saturated <- .top_covariate_matrix(saturated, "saturated", n)

  .top_validate_cutoff(subtype_cutoff)
  if (!is.null(start) && (!is.numeric(start) || any(!is.finite(start)))) {
    stop("`start` must be NULL or a finite numeric vector", call. = FALSE)
  }

  result <- TwoStageModel(
    y = y,
    baselineonly = baseline,
    additive = additive,
    pairwise.interaction = pairwise,
    saturated = saturated,
    missingTumorIndicator = missing_value,
    delta0 = start,
    cutoff = subtype_cutoff
  )
  if (!is.list(result)) {
    stop("TOP fitting failed without returning a model", call. = FALSE)
  }

  structure(
    result,
    class = c("top_model", "list"),
    call = match.call(),
    observations = n,
    tumor_markers = ncol(y) - 1L
  )
}

.top_validate_cutoff <- function(subtype_cutoff) {
  if (!is.numeric(subtype_cutoff) || length(subtype_cutoff) != 1L ||
      !is.finite(subtype_cutoff) || subtype_cutoff < 0 ||
      subtype_cutoff != floor(subtype_cutoff)) {
    stop("`subtype_cutoff` must be one non-negative integer", call. = FALSE)
  }
  invisible(subtype_cutoff)
}

#' @export
print.top_model <- function(x, ...) {
  cat("TOP two-stage polytomous model\n")
  cat("Observations:", attr(x, "observations"),
      " Tumor markers:", attr(x, "tumor_markers"), "\n")
  if (!is.null(x$loglikelihood)) {
    cat("Log-likelihood:", format(x$loglikelihood, digits = 6),
        " AIC:", format(x$AIC, digits = 6), "\n")
  }
  invisible(x)
}

#' @export
coef.top_model <- function(object, stage = c("second", "first"), ...) {
  stage <- match.arg(stage)
  if (stage == "second") object$delta else object$beta
}

#' @export
vcov.top_model <- function(object, stage = c("second", "first"), ...) {
  stage <- match.arg(stage)
  if (stage == "second") object$covariance.delta else object$covariance.beta
}

#' @export
summary.top_model <- function(object, ...) {
  structure(
    list(
      call = attr(object, "call"),
      second_stage = object$second.stage.test,
      global_tests = object$global.test,
      first_stage = object$first.stage.test,
      loglikelihood = object$loglikelihood,
      AIC = object$AIC
    ),
    class = "summary.top_model"
  )
}

#' @export
print.summary.top_model <- function(x, ...) {
  cat("TOP model summary\n")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("Log-likelihood:", format(x$loglikelihood, digits = 6),
      " AIC:", format(x$AIC, digits = 6), "\n")
  cat("\nGlobal tests:\n")
  print(x$global_tests)
  invisible(x)
}

#' Fit a TOP model with a custom second-stage design
#'
#' Use this interface when predefined baseline, additive, pairwise, or saturated
#' structures do not represent the desired subtype grouping.
#'
#' @param y Numeric outcome matrix as described in [top_fit()].
#' @param covariate Numeric vector or one-column matrix whose effect is mapped by
#'   `design`.
#' @param design Numeric matrix with one row per retained subtype and one column
#'   per custom second-stage effect.
#' @param adjustment Optional additive adjustment covariates.
#' @param missing_value Numeric missing-tumor code for cases, or `NULL`.
#' @param start Optional starting vector.
#' @param subtype_cutoff Subtypes with counts at or below this value are omitted.
#'
#' @return A named object of class `top_model`. Summary tables include both the
#'   custom covariate and any jointly fitted adjustment covariates.
#' @export
top_fit_custom <- function(y, covariate, design, adjustment = NULL,
                           missing_value = NULL, start = NULL,
                           subtype_cutoff = 10) {
  y <- .top_validate_outcome(y, missing_value)
  n <- nrow(y)
  covariate <- .top_covariate_matrix(covariate, "custom", n)
  if (ncol(covariate) != 1L) {
    stop("`covariate` must contain exactly one column", call. = FALSE)
  }
  adjustment <- .top_covariate_matrix(adjustment, "adjustment", n)
  .top_validate_cutoff(subtype_cutoff)
  if (!is.numeric(design) || length(dim(design)) != 2L ||
      ncol(design) == 0L || any(!is.finite(design))) {
    stop("`design` must be a finite numeric matrix with at least one column", call. = FALSE)
  }
  design <- as.matrix(design)
  expected_rows <- nrow(GenerateZstandard(
    y,
    missingTumorIndicator = missing_value,
    cutoff = subtype_cutoff
  ))
  if (nrow(design) != expected_rows) {
    stop("`design` must have ", expected_rows, " rows for the retained subtypes",
         call. = FALSE)
  }
  if (any(rowSums(abs(design)) == 0)) {
    stop("every retained subtype must map to at least one custom effect", call. = FALSE)
  }
  if (is.null(colnames(design))) {
    colnames(design) <- paste0("custom_effect", seq_len(ncol(design)))
  }

  fitting_function <- if (is.null(missing_value)) MvpolySelfDesign else EMmvpolySelfDesign
  result <- fitting_function(
    y = y,
    x.self.design = covariate,
    z.design = design,
    additive = adjustment,
    missingTumorIndicator = missing_value,
    delta0 = start,
    cutoff = subtype_cutoff
  )
  if (!is.list(result)) {
    stop("TOP fitting failed without returning a model", call. = FALSE)
  }
  result <- .top_complete_custom_summaries(
    result = result,
    y = y,
    covariate = covariate,
    design = design,
    adjustment = adjustment
  )
  structure(
    result,
    class = c("top_model", "list"),
    call = match.call(),
    observations = n,
    tumor_markers = ncol(y) - 1L,
    custom_design = design
  )
}

.top_complete_custom_summaries <- function(result, y, covariate, design,
                                           adjustment) {
  M <- nrow(result$z.standard)
  delta_index <- seq.int(M + 1L, length(result$delta))
  delta <- result$delta[delta_index]
  sigma <- result$covariance.delta[delta_index, delta_index, drop = FALSE]

  covariate_names <- colnames(covariate)
  block_sizes <- ncol(design)
  effect_names <- list(colnames(design))
  custom_blocks <- TRUE

  if (!is.null(adjustment)) {
    tumor_names <- colnames(y)[-1L]
    if (is.null(tumor_names)) {
      tumor_names <- as.character(seq_len(ncol(y) - 1L))
    }
    additive_names <- GenerateZDesignNamesAdditive(tumor_names)
    covariate_names <- c(covariate_names, colnames(adjustment))
    block_sizes <- c(block_sizes, rep(length(additive_names), ncol(adjustment)))
    effect_names <- c(effect_names, rep(list(additive_names), ncol(adjustment)))
    custom_blocks <- c(custom_blocks, rep(FALSE, ncol(adjustment)))
  }

  if (sum(block_sizes) != length(delta)) {
    stop("TOP returned an unexpected custom-model parameter layout", call. = FALSE)
  }

  block_end <- cumsum(block_sizes)
  block_start <- c(1L, block_end[-length(block_end)] + 1L)
  blocks <- Map(seq.int, block_start, block_end)
  second_stage_matrix <- matrix(
    NA_real_,
    nrow = max(block_sizes),
    ncol = length(block_sizes),
    dimnames = list(
      paste0("effect", seq_len(max(block_sizes))),
      covariate_names
    )
  )
  for (i in seq_along(blocks)) {
    second_stage_matrix[seq_len(block_sizes[i]), i] <- delta[blocks[[i]]]
  }
  attr(second_stage_matrix, "effect_names") <- effect_names

  standard_error <- sqrt(diag(sigma))
  result$second.stage.mat <- second_stage_matrix
  result$second.stage.test <- data.frame(
    Covariate = rep(covariate_names, block_sizes),
    SecondStageEffect = unlist(effect_names, use.names = FALSE),
    OddsRatio = round(exp(delta), 2L),
    `OddsRatio(95%CI low)` = round(exp(delta - 1.96 * standard_error), 2L),
    `OddsRatio(95%CI high)` = round(exp(delta + 1.96 * standard_error), 2L),
    Pvalue = IndividualHeterTest(delta, sigma),
    check.names = FALSE
  )

  result$global.test <- data.frame(
    Covariate = covariate_names,
    `Global test for association` = vapply(
      blocks,
      function(index) GlobalTestForAssoc(delta[index], sigma[index, index, drop = FALSE]),
      numeric(1)
    ),
    `Global test for heterogeneity` = vapply(
      seq_along(blocks),
      function(i) GlobalTestForHeter(
        delta[blocks[[i]]],
        sigma[blocks[[i]], blocks[[i]], drop = FALSE],
        self.design = custom_blocks[i]
      ),
      numeric(1)
    ),
    check.names = FALSE
  )

  subtype_names <- rownames(result$first.stage.mat)
  first_stage_matrix <- GenerateFirstStageMat(
    result$beta,
    covariate_names,
    subtype_names
  )
  first_stage_columns <- 1L + length(covariate_names)
  intercept_index <- 1L + first_stage_columns * (seq_len(M) - 1L)
  beta <- result$beta[-intercept_index]
  beta_sigma <- result$covariance.beta[
    -intercept_index,
    -intercept_index,
    drop = FALSE
  ]
  result$first.stage.mat <- first_stage_matrix
  result$first.stage.test <- FirstStageTest(
    beta,
    beta_sigma,
    M,
    first_stage_matrix
  )
  result
}
