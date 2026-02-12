#' Heatmap Visualization
#'
#' @name plot-heatmap
#' @keywords internal
NULL

#' Plot presence/absence heatmap
#'
#' Creates a heatmap visualization of the presence/absence matrix using
#' ComplexHeatmap. Rows are orthogroups, columns are assemblies.
#'
#' @param pa A `pa_matrix` object.
#' @param cluster_rows Logical. Whether to cluster rows (default: TRUE).
#' @param cluster_cols Logical. Whether to cluster columns (default: TRUE).
#' @param show_row_names Logical. Whether to show row names (default: FALSE).
#' @param show_col_names Logical. Whether to show column names (default: TRUE).
#' @param color Character vector or color function. Colors for the heatmap.
#'   If NULL, auto-detected based on pa$type.
#' @param row_annotation A data.frame with row names matching orthogroup IDs.
#'   Columns become annotation tracks on the left side.
#' @param col_annotation A data.frame with row names matching assembly names.
#'   Columns become annotation tracks on the top.
#' @param distance_method Character. Distance metric for clustering:
#'   "binary" (default), "jaccard", or "bray-curtis".
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap().
#'
#' @return A ComplexHeatmap Heatmap object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_heatmap(pa)
#' plot_heatmap(pa, cluster_rows = FALSE)
#' }
plot_heatmap <- function(pa,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_row_names = FALSE,
                         show_col_names = TRUE,
                         color = NULL,
                         row_annotation = NULL,
                         col_annotation = NULL,
                         distance_method = "binary",
                         ...) {

  # Validate input
  if (!inherits(pa, "pa_matrix")) {
    cli::cli_abort("{.arg pa} must be a {.cls pa_matrix} object")
  }

  # Validate distance_method
  valid_methods <- c("binary", "jaccard", "bray-curtis")
  if (!distance_method %in% valid_methods) {
    cli::cli_abort(
      "{.arg distance_method} must be one of {.val {valid_methods}}, not {.val {distance_method}}"
    )
  }

  # Process row annotation
  row_ha <- NULL
  if (!is.null(row_annotation)) {
    row_ann_processed <- validate_annotation(
      row_annotation,
      rownames(pa$matrix),
      "row"
    )
    if (!is.null(row_ann_processed) && nrow(row_ann_processed) > 0) {
      row_ha <- ComplexHeatmap::rowAnnotation(df = row_ann_processed)
    }
  }

  # Process column annotation
  col_ha <- NULL
  if (!is.null(col_annotation)) {
    col_ann_processed <- validate_annotation(
      col_annotation,
      colnames(pa$matrix),
      "column"
    )
    if (!is.null(col_ann_processed) && nrow(col_ann_processed) > 0) {
      col_ha <- ComplexHeatmap::HeatmapAnnotation(df = col_ann_processed)
    }
  }

  # Auto-detect color scheme if not provided
  if (is.null(color)) {
    max_val <- max(pa$matrix, na.rm = TRUE)
    min_val <- min(pa$matrix, na.rm = TRUE)

    color <- switch(pa$type,
      "binary" = c("0" = "white", "1" = "steelblue"),
      "count" = {
        # Handle case where all values are the same
        if (max_val == min_val) {
          c("white", "steelblue")
        } else {
          circlize::colorRamp2(
            c(0, max_val),
            c("white", "steelblue")
          )
        }
      },
      "score" = {
        # Handle case where all values are the same
        if (max_val == min_val || is.na(max_val)) {
          c("white", "firebrick")
        } else {
          circlize::colorRamp2(
            c(0, max_val),
            c("white", "firebrick")
          )
        }
      }
    )
  }

  # Handle edge cases: disable clustering for single row/column
  n_rows <- nrow(pa$matrix)
  n_cols <- ncol(pa$matrix)

  if (n_rows < 2) {
    cluster_rows <- FALSE
  }
  if (n_cols < 2) {
    cluster_cols <- FALSE
  }

  # Pre-compute clustering with specified distance method
  # Handle NA values (score matrices have NAs for absent orthogroups)
  m_clean <- pa$matrix
  if (any(is.na(m_clean))) {
    m_clean[is.na(m_clean)] <- 0
  }

  if (cluster_rows && n_rows >= 2) {
    row_dist <- calculate_distance(m_clean, distance_method)
    cluster_rows <- stats::hclust(row_dist)
  }
  if (cluster_cols && n_cols >= 2) {
    col_dist <- calculate_distance(t(m_clean), distance_method)
    cluster_cols <- stats::hclust(col_dist)
  }

  # Create heatmap
  ComplexHeatmap::Heatmap(
    pa$matrix,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_row_names,
    show_column_names = show_col_names,
    col = color,
    name = pa$type,
    na_col = "grey90",
    left_annotation = row_ha,
    top_annotation = col_ha,
    ...
  )
}

#' Calculate distance matrix using specified method
#'
#' @param mat Numeric matrix (rows are observations)
#' @param method One of "binary", "jaccard", "bray-curtis"
#' @return A dist object
#' @keywords internal
calculate_distance <- function(mat, method) {
  switch(method,
    "binary" = stats::dist(mat, method = "binary"),
    "jaccard" = jaccard_dist(mat),
    "bray-curtis" = bray_curtis_dist(mat)
  )
}

#' Calculate Jaccard distance
#'
#' Jaccard dissimilarity: 1 - |A ∩ B| / |A ∪ B|
#' For binary data: 1 - (number of shared 1s) / (number where either is 1)
#'
#' @param mat Numeric matrix (rows are observations)
#' @return A dist object
#' @keywords internal
jaccard_dist <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      a <- mat[i, ]
      b <- mat[j, ]

      # For binary: intersection = both > 0, union = either > 0
      intersection <- sum(a > 0 & b > 0)
      union <- sum(a > 0 | b > 0)

      if (union == 0) {
        d[i, j] <- 0
      } else {
        d[i, j] <- 1 - (intersection / union)
      }
      d[j, i] <- d[i, j]
    }
  }

  stats::as.dist(d)
}

#' Calculate Bray-Curtis distance
#'
#' Bray-Curtis dissimilarity: sum(|xi - yi|) / sum(xi + yi)
#'
#' @param mat Numeric matrix (rows are observations)
#' @return A dist object
#' @keywords internal
bray_curtis_dist <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      a <- mat[i, ]
      b <- mat[j, ]

      numerator <- sum(abs(a - b))
      denominator <- sum(a + b)

      if (denominator == 0) {
        d[i, j] <- 0
      } else {
        d[i, j] <- numerator / denominator
      }
      d[j, i] <- d[i, j]
    }
  }

  stats::as.dist(d)
}

#' Validate and process annotation data frame
#'
#' @param annotation Data frame with row names
#' @param expected_names Character vector of expected row names
#' @param type "row" or "column" for error messages
#' @return Processed annotation data frame (reordered to match expected_names, with NA for missing)
#' @keywords internal
validate_annotation <- function(annotation, expected_names, type) {
  if (!is.data.frame(annotation)) {
    cli::cli_abort("{.arg {type}_annotation} must be a data.frame")
  }

  ann_names <- rownames(annotation)
  matching <- intersect(expected_names, ann_names)

  if (length(matching) == 0) {
    cli::cli_warn(
      "No {type} annotation names match the matrix {type}s. Annotation will be ignored."
    )
    return(NULL)
  }

  if (length(matching) < length(expected_names)) {
    missing <- setdiff(expected_names, ann_names)
    cli::cli_warn(
      "{length(missing)} {type}(s) have no matching annotation: {.val {head(missing, 3)}}{if (length(missing) > 3) '...' else ''}"
    )
  }

  # Build result data frame with expected_names as rows (missing = NA)
  result <- data.frame(row.names = expected_names)
  for (col_name in names(annotation)) {
    col_values <- annotation[[col_name]]
    new_col <- rep(NA, length(expected_names))

    # Match by row names
    idx <- match(expected_names, ann_names)
    has_match <- !is.na(idx)
    new_col[has_match] <- col_values[idx[has_match]]

    # Preserve original type
    if (is.factor(col_values)) {
      new_col <- factor(new_col, levels = levels(col_values))
    } else if (is.character(col_values)) {
      new_col <- as.character(new_col)
    } else if (is.numeric(col_values)) {
      new_col <- as.numeric(new_col)
      # Handle edge case: if only one unique non-NA value, convert to character
      # to avoid colorRamp2 "need at least two distinct break values" error
      unique_vals <- unique(new_col[!is.na(new_col)])
      if (length(unique_vals) <= 1) {
        new_col <- as.character(new_col)
      }
    }

    result[[col_name]] <- new_col
  }

  result
}
