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
#' @section Styles:
#' `style` selects a ready-made look. All three use the same layout rules
#' (light grid between cells on matrices of up to 200 rows, slim dendrograms,
#' angled column labels, legends labelled "Absent"/"Present") and differ in
#' colour and in what the rows say:
#'
#' * `"clean"` (default): one hue on a light grey background. Binary and count
#'   matrices are blue; score matrices are orange so that they are never
#'   mistaken for presence/absence.
#' * `"pan"`: rows are split into **Core** (present in every assembly),
#'   **Accessory** (present in more than one but not all) and **Unique**
#'   (present in exactly one) blocks, each with a coloured side bar. On binary
#'   matrices the cells themselves take the class colour; on count and score
#'   matrices the cells keep the magnitude ramp and the class is shown by the
#'   split. Note that `plot_pan_structure()` further divides Accessory into
#'   Accessory and Rare using thresholds; the heatmap does not.
#' * `"mono"`: black on white with grey grid lines, for print and greyscale
#'   figures.
#'
#' `color` overrides the cell colours of any style. Count matrices whose
#' maximum is 1 are drawn with the binary colours.
#'
#' @section Customising further:
#' The return value is a `ComplexHeatmap::Heatmap` object, so anything that
#' package can do applies. Arguments in `...` are passed to
#' [ComplexHeatmap::Heatmap()] and override the style's settings, e.g.
#' `column_title = "Effectors"`, `rect_gp = grid::gpar(col = NA)` to remove
#' the grid, or `column_names_rot = 90`. See `vignette("styling-heatmaps")`.
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
#' @param style Character. One of `"clean"` (default), `"pan"` or `"mono"`.
#'   See "Styles".
#' @param cluster_rows Logical. Whether to cluster rows (default: TRUE).
#' @param cluster_cols Logical. Whether to cluster columns (default: TRUE).
#' @param show_row_names Logical. Whether to show row names (default: FALSE).
#' @param show_col_names Logical. Whether to show column names (default: TRUE).
#' @param color Character vector or color function. Colors for the heatmap.
#'   If NULL, chosen from `style` and `pa$type`.
#' @param row_annotation A data.frame with row names matching orthogroup IDs.
#'   Columns become annotation tracks on the left side.
#' @param col_annotation A data.frame with row names matching assembly names.
#'   Columns become annotation tracks on the top.
#' @param distance_method Character. Distance metric for clustering:
#'   "binary" (default), "jaccard", or "bray-curtis".
#' @param max_cluster_rows Integer. Above this many rows, row clustering is
#'   done on distinct patterns rather than individual rows (default 5000).
#'   See "Large matrices".
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap(). These
#'   take precedence over the style's settings.
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
#' # Core / accessory / unique blocks
#' plot_heatmap(pa, style = "pan")
#'
#' # Greyscale for print, no row clustering
#' plot_heatmap(pa, style = "mono", cluster_rows = FALSE, distance_method = "jaccard")
plot_heatmap <- function(pa,
                         style = c("clean", "pan", "mono"),
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

  style <- rlang::arg_match(style)
  st <- heatmap_style(style)

  # Validate distance_method
  valid_methods <- c("binary", "jaccard", "bray-curtis")
  if (!distance_method %in% valid_methods) {
    cli::cli_abort(
      "{.arg distance_method} must be one of {.val {valid_methods}}, not {.val {distance_method}}"
    )
  }

  mat <- pa$matrix
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)

  # NA-free copy for distances and classification (score matrices have NAs)
  m_clean <- mat
  if (any(is.na(m_clean))) {
    m_clean[is.na(m_clean)] <- 0
  }

  # Orthogroup classes for the pan style
  og_class <- NULL
  if (style == "pan") {
    og_class <- classify_rows(m_clean)
  }

  # Process row annotation
  row_ann_processed <- NULL
  if (!is.null(row_annotation)) {
    row_ann_processed <- validate_annotation(row_annotation, rownames(mat), "row")
    if (!is.null(row_ann_processed) && nrow(row_ann_processed) == 0) {
      row_ann_processed <- NULL
    }
  }
  row_ha <- NULL
  if (style == "pan") {
    block_cols <- st$class_colors[levels(og_class)]
    ha_args <- list(
      class = ComplexHeatmap::anno_block(
        gp = grid::gpar(fill = block_cols, col = NA),
        width = grid::unit(2, "mm"),
        which = "row"
      ),
      show_annotation_name = FALSE
    )
    if (!is.null(row_ann_processed)) {
      ha_args$df <- row_ann_processed
      ha_args$show_annotation_name <- c(FALSE, rep(TRUE, ncol(row_ann_processed)))
    }
    row_ha <- do.call(ComplexHeatmap::rowAnnotation, ha_args)
  } else if (!is.null(row_ann_processed)) {
    row_ha <- ComplexHeatmap::rowAnnotation(df = row_ann_processed)
  }

  # Process column annotation
  col_ha <- NULL
  if (!is.null(col_annotation)) {
    col_ann_processed <- validate_annotation(col_annotation, colnames(mat), "column")
    if (!is.null(col_ann_processed) && nrow(col_ann_processed) > 0) {
      col_ha <- ComplexHeatmap::HeatmapAnnotation(df = col_ann_processed)
    }
  }

  # Colours and legend
  legend_param <- list(
    title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9),
    border = st$legend_border
  )
  legend_name <- switch(pa$type, binary = "Orthogroup", count = "Copies", score = "Score")
  max_val <- suppressWarnings(max(mat, na.rm = TRUE))
  min_val <- suppressWarnings(min(mat, na.rm = TRUE))
  draw_mat <- mat

  if (is.null(color)) {
    binary_like <- pa$type == "binary" ||
      (pa$type == "count" && is.finite(max_val) && max_val <= 1)

    if (binary_like && style == "pan") {
      # Cells coloured by orthogroup class: 0 absent, then class index
      draw_mat <- (m_clean > 0) * as.integer(og_class)
      storage.mode(draw_mat) <- "integer"
      used <- levels(og_class)
      color <- c("0" = st$absent, stats::setNames(st$class_colors[used], seq_along(used)))
      legend_param$at <- c(seq_along(used), 0)
      legend_param$labels <- c(used, "Absent")
      legend_name <- "Orthogroup"
    } else if (binary_like) {
      color <- c("0" = st$absent, "1" = st$present)
      legend_param$at <- c(0, 1)
      legend_param$labels <- c("Absent", "Present")
      legend_name <- "Orthogroup"
    } else if (pa$type == "count") {
      color <- circlize::colorRamp2(c(0, 1, max_val), st$count_ramp)
    } else {
      # score: absent cells are NA and drawn in the absent colour
      if (!is.finite(max_val) || max_val == min_val) {
        color <- c("0" = st$absent, "1" = st$score_ramp[2])
        draw_mat <- (m_clean > 0) * 1L
        legend_param$at <- c(0, 1)
        legend_param$labels <- c("Absent", "Present")
      } else {
        color <- circlize::colorRamp2(c(min_val, max_val), st$score_ramp)
      }
    }
  }

  # Handle edge cases: disable clustering for single row/column
  if (n_rows < 2) {
    cluster_rows <- FALSE
  }
  if (n_cols < 2) {
    cluster_cols <- FALSE
  }

  # Row clustering. With a row split (pan style) ComplexHeatmap must cluster
  # within each slice itself, so we hand it the distance function instead of
  # a precomputed hclust.
  row_order <- NULL
  clustering_distance_rows <- NULL
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
    } else if (style == "pan") {
      cluster_rows <- TRUE
      clustering_distance_rows <- function(m) {
        m[is.na(m)] <- 0
        calculate_distance(m, distance_method)
      }
    } else {
      row_dist <- calculate_distance(m_clean, distance_method)
      cluster_rows <- stats::hclust(row_dist)
    }
  }
  if (isTRUE(cluster_cols) && n_cols >= 2) {
    col_dist <- calculate_distance(t(m_clean), distance_method)
    cluster_cols <- stats::hclust(col_dist)
  }

  # Assemble arguments: style defaults first, then user `...` overrides
  ht_args <- list(
    draw_mat,
    name = legend_name,
    col = color,
    na_col = st$absent,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_row_names,
    show_column_names = show_col_names,
    left_annotation = row_ha,
    top_annotation = col_ha,
    heatmap_legend_param = legend_param,
    rect_gp = if (n_rows <= 200) st$grid_gp else grid::gpar(col = NA),
    border = st$border,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 9, col = st$label_col),
    row_names_gp = grid::gpar(fontsize = 8, col = st$label_col),
    row_dend_width = grid::unit(10, "mm"),
    column_dend_height = grid::unit(8, "mm"),
    row_dend_gp = st$dend_gp,
    column_dend_gp = st$dend_gp
  )
  if (!is.null(row_order)) {
    ht_args$row_order <- row_order
  }
  if (!is.null(clustering_distance_rows)) {
    ht_args$clustering_distance_rows <- clustering_distance_rows
  }
  if (style == "pan") {
    ht_args$row_split <- og_class
    ht_args$cluster_row_slices <- FALSE
    ht_args$row_title_rot <- 0
    ht_args$row_title_gp <- grid::gpar(fontsize = 9, fontface = "bold", col = st$label_col)
    ht_args$row_gap <- grid::unit(1.5, "mm")
  }

  dots <- list(...)
  if (length(dots)) {
    ht_args[names(dots)] <- dots
  }
  do.call(ComplexHeatmap::Heatmap, ht_args)
}

#' Heatmap style parameters
#'
#' Colours and graphical parameters for each `plot_heatmap()` style. Hues are
#' taken from a palette validated for colour-vision deficiency: Core blue,
#' Accessory green and Unique orange are mutually distinguishable under the
#' common deficiencies.
#'
#' @param style One of "clean", "pan", "mono".
#' @return A named list.
#' @keywords internal
heatmap_style <- function(style) {
  absent_light <- "#f0efec"
  ink <- "#1a1a19"
  switch(style,
    clean = ,
    pan = list(
      absent = absent_light,
      present = "#2a78d6",
      count_ramp = c(absent_light, "#9ec5f4", "#0d366b"),
      score_ramp = c("#fde2d3", "#b3400f"),
      class_colors = c(Core = "#2a78d6", Accessory = "#1baf7a",
                       Unique = "#eb6834", Absent = "#c3c2b7"),
      grid_gp = grid::gpar(col = "white", lwd = 0.8),
      dend_gp = grid::gpar(col = "#8a8985", lwd = 0.6),
      label_col = "#52514e",
      border = FALSE,
      legend_border = NA
    ),
    mono = list(
      absent = "white",
      present = ink,
      count_ramp = c("white", "#c3c2b7", ink),
      score_ramp = c("#d4d3cf", ink),
      class_colors = c(Core = ink, Accessory = "#6f6e6a",
                       Unique = "#b5b4ae", Absent = "#e5e4e0"),
      grid_gp = grid::gpar(col = "#c3c2b7", lwd = 0.5),
      dend_gp = grid::gpar(col = ink, lwd = 0.6),
      label_col = "black",
      border = ink,
      legend_border = ink
    )
  )
}

#' Classify orthogroups by how many assemblies they occur in
#'
#' Core: present in every column; Unique: present in exactly one; Accessory:
#' anything in between. Rows with no presence at all (possible after
#' subsetting) are labelled Absent. Unused levels are dropped so that the
#' heatmap split only shows classes that occur.
#'
#' @param mat Numeric matrix without NAs.
#' @return A factor with levels in the order Core, Accessory, Unique, Absent.
#' @keywords internal
classify_rows <- function(mat) {
  n <- rowSums(mat > 0)
  k <- ncol(mat)
  cls <- ifelse(n == 0, "Absent",
         ifelse(n == k, "Core",
         ifelse(n == 1, "Unique", "Accessory")))
  droplevels(factor(cls, levels = c("Core", "Accessory", "Unique", "Absent")))
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
