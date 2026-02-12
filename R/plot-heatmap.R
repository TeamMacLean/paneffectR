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
                         ...) {
  # Validate input
 if (!inherits(pa, "pa_matrix")) {
    cli::cli_abort("{.arg pa} must be a {.cls pa_matrix} object")
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

  # Handle NA values for clustering (score matrices have NAs for absent orthogroups)
  # Pre-compute clustering with NA-safe distance
  has_na <- any(is.na(pa$matrix))
  if (has_na && (cluster_rows || cluster_cols)) {
    # Replace NA with 0 for distance calculation
    m_clean <- pa$matrix
    m_clean[is.na(m_clean)] <- 0

    if (cluster_rows && n_rows >= 2) {
      row_dist <- stats::dist(m_clean, method = "binary")
      cluster_rows <- stats::hclust(row_dist)
    }
    if (cluster_cols && n_cols >= 2) {
      col_dist <- stats::dist(t(m_clean), method = "binary")
      cluster_cols <- stats::hclust(col_dist)
    }
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
    ...
  )
}
