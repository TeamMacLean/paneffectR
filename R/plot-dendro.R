#' Dendrogram Visualization
#'
#' @name plot-dendro
#' @keywords internal
NULL

#' Plot assembly clustering dendrogram
#'
#' Clusters assemblies by shared orthogroup content and displays as a dendrogram.
#'
#' @param pa A `pa_matrix` object.
#' @param distance_method Character. Distance metric: "jaccard" (default),
#'   "binary", or "bray-curtis".
#' @param cluster_method Character. Clustering method for hclust
#'   (default: "complete"). Options: "complete", "single", "average", "ward.D2".
#' @param labels Logical. Whether to show assembly labels (default: TRUE).
#' @param hang Numeric. Fraction of plot height to hang labels below leaves
#'   (default: 0.1).
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_dendro(pa)
#' plot_dendro(pa, distance_method = "bray-curtis")
#' plot_dendro(pa, cluster_method = "ward.D2")
#' }
plot_dendro <- function(pa,
                        distance_method = "jaccard",
                        cluster_method = "complete",
                        labels = TRUE,
                        hang = 0.1) {
  # Validate input
  if (!inherits(pa, "pa_matrix")) {
    cli::cli_abort("{.arg pa} must be a {.cls pa_matrix} object")
  }

  # Validate distance_method
  valid_methods <- c("jaccard", "binary", "bray-curtis")
  if (!distance_method %in% valid_methods) {
    cli::cli_abort(
      "{.arg distance_method} must be one of {.val {valid_methods}}, not {.val {distance_method}}"
    )
  }

  # Need at least 2 assemblies to cluster
  if (ncol(pa$matrix) < 2) {
    cli::cli_abort("Need at least 2 assemblies to create dendrogram")
  }

  # Calculate distance between assemblies (transpose so assemblies are rows)
  mat_t <- t(pa$matrix)
  mat_t[is.na(mat_t)] <- 0

  d <- calculate_assembly_distance(mat_t, distance_method)

  # Hierarchical clustering
  hc <- stats::hclust(d, method = cluster_method)

  # Convert to dendrogram data for ggplot2
  dend <- stats::as.dendrogram(hc)
  dend_data <- dendro_to_data(dend)

  # Build ggplot
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dend_data$segments,
      ggplot2::aes(
        x = .data$x, y = .data$y,
        xend = .data$xend, yend = .data$yend
      ),
      color = "black"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Distance",
      title = paste0("Assembly clustering (", distance_method, ")")
    )

  # Add labels if requested
  if (labels && !is.null(dend_data$labels) && nrow(dend_data$labels) > 0) {
    y_offset <- max(dend_data$segments$y) * hang
    p <- p + ggplot2::geom_text(
      data = dend_data$labels,
      ggplot2::aes(x = .data$x, y = -y_offset, label = .data$label),
      hjust = 0.5, vjust = 1, size = 3, angle = 45
    )
  }

  p
}

#' Calculate distance between assemblies
#'
#' @param mat Matrix with assemblies as rows, orthogroups as columns
#' @param method Distance method
#' @return dist object
#' @keywords internal
calculate_assembly_distance <- function(mat, method) {
  switch(method,
    "binary" = stats::dist(mat, method = "binary"),
    "jaccard" = jaccard_dist(mat),
    "bray-curtis" = bray_curtis_dist(mat)
  )
}

#' Convert dendrogram to data frames for ggplot2
#'
#' @param dend A dendrogram object
#' @return List with segments and labels data frames
#' @keywords internal
dendro_to_data <- function(dend) {
  # Get segment data by traversing the dendrogram
  segments <- list()
  labels_data <- list()

  # Helper function to extract segments recursively
  extract_segments <- function(d, x_center = 0) {
    if (is.leaf(d)) {
      # Leaf node - collect label
      labels_data[[length(labels_data) + 1]] <<- data.frame(
        x = x_center,
        y = 0,
        label = attr(d, "label"),
        stringsAsFactors = FALSE
      )
      return(x_center)
    }

    # Get positions of children
    members <- attr(d, "members")
    height <- attr(d, "height")

    # Process left and right children
    left <- d[[1]]
    right <- d[[2]]

    left_members <- attr(left, "members")
    right_members <- attr(right, "members")

    # Calculate x positions for children
    left_x <- x_center - right_members / 2
    right_x <- x_center + left_members / 2

    # Process children
    left_x_actual <- extract_segments(left, left_x)
    right_x_actual <- extract_segments(right, right_x)

    # Get heights of children
    left_height <- if (is.leaf(left)) 0 else attr(left, "height")
    right_height <- if (is.leaf(right)) 0 else attr(right, "height")

    # Add segments: vertical from children to parent height, horizontal at parent
    # Left vertical
    segments[[length(segments) + 1]] <<- data.frame(
      x = left_x_actual, y = left_height,
      xend = left_x_actual, yend = height
    )
    # Right vertical
    segments[[length(segments) + 1]] <<- data.frame(
      x = right_x_actual, y = right_height,
      xend = right_x_actual, yend = height
    )
    # Horizontal at parent level
    segments[[length(segments) + 1]] <<- data.frame(
      x = left_x_actual, y = height,
      xend = right_x_actual, yend = height
    )

    return(x_center)
  }

  extract_segments(dend, x_center = attr(dend, "members") / 2)

  list(
    segments = do.call(rbind, segments),
    labels = if (length(labels_data) > 0) do.call(rbind, labels_data) else NULL
  )
}
