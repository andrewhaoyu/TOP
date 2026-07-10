#' Fit detailed TOP models for selected follow-up variants
#'
#' Fits one TOP model per selected variant and returns variant-specific global
#' tests and second-stage effect estimates. This is intended for the limited
#' set of variants selected after a score scan, not for genome-wide fitting.
#'
#' @param y Numeric outcome matrix as described in [top_fit()].
#' @param variants Numeric vector, matrix, or data frame with selected variants
#'   in columns.
#' @param adjustment Optional numeric adjustment covariates. These use the
#'   additive second-stage structure.
#' @param variant_structure Second-stage structure for each tested variant.
#' @param missing_value Numeric missing-tumor code for cases, or `NULL`.
#' @param start Optional starting vector passed to [top_fit()].
#' @param subtype_cutoff Subtypes with counts at or below this value are omitted.
#' @param on_error Whether a failed fit should be recorded or stop the run.
#' @param keep_models Keep the fitted `top_model` objects in the result.
#'
#' @return A `top_followup` object with global tests, effect estimates, optional
#'   fitted models, and recorded errors.
#' @export
top_followup <- function(y, variants, adjustment = NULL,
                         variant_structure = c("additive", "baseline", "pairwise", "saturated"),
                         missing_value = NULL, start = NULL, subtype_cutoff = 10,
                         on_error = c("record", "stop"), keep_models = TRUE) {
  variant_structure <- match.arg(variant_structure)
  on_error <- match.arg(on_error)
  y <- .top_validate_outcome(y, missing_value)
  variants <- .top_variant_matrix(variants, nrow(y))
  adjustment <- .top_covariate_matrix(adjustment, "adjustment", nrow(y))
  .top_validate_cutoff(subtype_cutoff)

  models <- vector("list", ncol(variants))
  names(models) <- colnames(variants)
  result_rows <- vector("list", ncol(variants))
  effect_rows <- vector("list", ncol(variants))
  error_rows <- vector("list", ncol(variants))

  for (index in seq_len(ncol(variants))) {
    variant_name <- colnames(variants)[index]
    variant <- variants[, index, drop = FALSE]
    fit_arguments <- list(
      y = y,
      additive = adjustment,
      missing_value = missing_value,
      start = start,
      subtype_cutoff = subtype_cutoff
    )
    if (variant_structure == "additive") {
      fit_arguments$additive <- cbind(variant, adjustment)
    } else {
      fit_arguments[[variant_structure]] <- variant
    }

    fit <- tryCatch(do.call(top_fit, fit_arguments), error = function(error) error)
    if (inherits(fit, "error")) {
      if (on_error == "stop") stop(fit)
      message_text <- conditionMessage(fit)
      result_rows[[index]] <- data.frame(
        variant = variant_name,
        global_association_p = NA_real_,
        global_heterogeneity_p = NA_real_,
        status = message_text,
        stringsAsFactors = FALSE
      )
      error_rows[[index]] <- data.frame(
        variant = variant_name, error = message_text, stringsAsFactors = FALSE
      )
      next
    }

    models[[index]] <- fit
    global <- fit$global.test[fit$global.test$Covariate == variant_name, , drop = FALSE]
    effects <- fit$second.stage.test[
      fit$second.stage.test$Covariate == variant_name, , drop = FALSE
    ]
    if (nrow(global) != 1L) {
      stop("TOP could not identify variant-specific global tests for `", variant_name, "`",
           call. = FALSE)
    }
    result_rows[[index]] <- data.frame(
      variant = variant_name,
      global_association_p = global[["Global test for association"]],
      global_heterogeneity_p = global[["Global test for heterogeneity"]],
      status = "ok",
      stringsAsFactors = FALSE
    )
    effect_rows[[index]] <- effects
  }

  results <- do.call(rbind, result_rows)
  effects <- if (all(vapply(effect_rows, is.null, logical(1)))) {
    data.frame()
  } else {
    do.call(rbind, effect_rows[!vapply(effect_rows, is.null, logical(1))])
  }
  errors <- if (all(vapply(error_rows, is.null, logical(1)))) {
    data.frame(variant = character(), error = character())
  } else {
    do.call(rbind, error_rows[!vapply(error_rows, is.null, logical(1))])
  }

  structure(
    list(
      results = results,
      effects = effects,
      models = if (isTRUE(keep_models)) models else NULL,
      errors = errors,
      variant_structure = variant_structure
    ),
    class = "top_followup"
  )
}

#' @export
print.top_followup <- function(x, ...) {
  cat("TOP follow-up models\n")
  cat("Variant structure:", x$variant_structure,
      " Variants:", nrow(x$results), "\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.top_followup <- function(x, ...) x$results
