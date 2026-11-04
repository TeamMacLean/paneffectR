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
#' @section Large matrices:
#' Clustering every row of a matrix with tens of thousands of orthogroups is
#' slow and the resulting dendrogram is unreadable. When the matrix has more
#' than `max_cluster_rows` rows and `cluster_rows = TRUE`, rows are instead
#' ordered by clustering the distinct presence/absence patterns (at most
#' `2^n_assemblies`), which gives the same grouping in a fraction of the time
#' but draws no row dendrogram. For a readable plot, reduce the matrix first,
#' e.g. with `build_pa_matrix(..., exclude_singletons = TRUE)` or a
#' `score_threshold`.
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
#' @param max_cluster_rows Integer. Above this many rows, row clustering is
#'   done on distinct patterns rather than individual rows (default 5000).
#'   See "Large matrices".
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap().
#'
#' @return A ComplexHeatmap Heatmap object.
#' @export
#'
#' @examples
#' visual_dir <- system.file("testdata", "visual", package = "paneffectR")
#' pa <- readRDS(file.path(visual_dir, "pa_binary.rds"))
#' ht <- plot_heatmap(pa)
#' ComplexHeatmap::draw(ht)
#'
#' plot_heatmap(pa, cluster_rows = FALSE, distance_method = "jaccard")
plot_heatmap <- function(pa,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_row_names = FALSE,
                         show_col_names = TRUE,
                         color = NULL,
                         row_annotation = NULL,
                         col_annotation = NULL,
                         distance_method = "binary",
                         max_cluster_rows = 5000,
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

  row_order <- NULL
  if (isTRUE(cluster_rows) && n_rows >= 2) {
    if (n_rows > max_cluster_rows) {
      # Too many rows for a full dendrogram: order by clustering distinct patterns
      row_order <- pattern_row_order(m_clean, distance_method)
      n_patterns <- length(unique(row_order$pattern))
      cli::cli_inform(c(
        "i" = "{n_rows} rows exceed {.arg max_cluster_rows} ({max_cluster_rows}): rows are ordered by clustering their {n_patterns} distinct presence/absence pattern{?s} and no row dendrogram is drawn.",
        "i" = "Reduce rows with {.code exclude_singletons = TRUE} or a {.arg score_threshold} in {.fn build_pa_matrix}, or raise {.arg max_cluster_rows}."
      ))
      row_order <- row_order$order
      cluster_rows <- FALSE
    } else {
      row_dist <- calculate_distance(m_clean, distance_method)
      cluster_rows <- stats::hclust(row_dist)
    }
  }
  if (isTRUE(cluster_cols) && n_cols >= 2) {
    col_dist <- calculate_distance(t(m_clean), distance_method)
    cluster_cols <- stats::hclust(col_dist)
  }

  # Create heatmap
  ht_args <- list(
    pa$matrix,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_row_names,
    show_column_names = show_col_names,
    col = color,
    name = pa$type,
    na_col = "grey90",
    left_annotation = row_ha,
    top_annotation = col_ha
  )
  if (!is.null(row_order)) {
    ht_args$row_order <- row_order
  }
  do.call(ComplexHeatmap::Heatmap, c(ht_args, list(...)))
}

#' Order rows by clustering their distinct presence/absence patterns
#'
#' Rows with an identical pattern across columns are grouped; the distinct
#' patterns are clustered with `calculate_distance()` + `hclust()`, and every
#' row is placed according to its pattern's position in that clustering.
#'
#' @param mat Numeric matrix without NAs (rows are observations).
#' @param method Distance method passed to `calculate_distance()`.
#' @return A list with `order` (integer row order) and `pattern` (integer
#'   pattern index per row, in original row order).
#' @keywords internal
pattern_row_order <- function(mat, method) {
  key <- do.call(paste, c(as.data.frame(mat), sep = ","))
  pattern <- match(key, unique(key))
  patterns <- mat[!duplicated(key), , drop = FALSE]
  rownames(patterns) <- NULL

  if (nrow(patterns) >= 2) {
    pattern_order <- stats::hclust(calculate_distance(patterns, method))$order
  } else {
    pattern_order <- seq_len(nrow(patterns))
  }

  list(
    order = order(match(pattern, pattern_order)),
    pattern = pattern
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
  b <- (mat > 0) * 1
  intersection <- b %*% t(b)
  rs <- rowSums(b)
  union <- outer(rs, rs, "+") - intersection

  d <- 1 - intersection / union
  d[union == 0] <- 0
  dimnames(d) <- NULL
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
  numerator <- as.matrix(stats::dist(mat, method = "manhattan"))
  rs <- rowSums(mat)
  denominator <- outer(rs, rs, "+")

  d <- numerator / denominator
  d[denominator == 0] <- 0
  dimnames(d) <- NULL
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
