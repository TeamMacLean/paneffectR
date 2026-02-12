#' Score Distribution Visualization
#'
#' @name plot-scores
#' @keywords internal
NULL

#' Plot effector score distributions
#'
#' Creates density plots or histograms of effector scores from a protein collection.
#'
#' @param proteins A `protein_collection` object with score data.
#' @param score_column Character. Column name for scores (default: "custom_score").
#' @param by_assembly Logical. Whether to facet by assembly (default: FALSE).
#'   If TRUE, creates separate panels for each assembly.
#' @param threshold Numeric or NULL. If specified, draws a vertical line at this
#'   value to show where a score cutoff would fall.
#' @param plot_type Character. Plot type: "density" (default) or "histogram".
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_scores(proteins)
#' plot_scores(proteins, threshold = 5.0)
#' plot_scores(proteins, by_assembly = TRUE)
#' plot_scores(proteins, plot_type = "histogram")
#' }
plot_scores <- function(proteins,
                        score_column = "custom_score",
                        by_assembly = FALSE,
                        threshold = NULL,
                        plot_type = "density") {
  # Validate input
  if (!inherits(proteins, "protein_collection")) {
    cli::cli_abort("{.arg proteins} must be a {.cls protein_collection} object")
  }

  # Validate plot_type
  valid_types <- c("density", "histogram")
  if (!plot_type %in% valid_types) {
    cli::cli_abort(
      "{.arg plot_type} must be one of {.val {valid_types}}, not {.val {plot_type}}"
    )
  }

  # Extract scores from all assemblies
  score_data <- extract_scores(proteins, score_column)

  if (is.null(score_data) || nrow(score_data) == 0) {
    cli::cli_abort("No scores found in column {.val {score_column}}")
  }

  # Remove NAs
  score_data <- score_data[!is.na(score_data$score), ]

  if (nrow(score_data) == 0) {
    cli::cli_abort("All scores are NA in column {.val {score_column}}")
  }

  # Build base plot
  p <- ggplot2::ggplot(score_data, ggplot2::aes(x = .data$score))

  # Add geom based on plot_type
  if (plot_type == "density") {
    if (by_assembly) {
      p <- p + ggplot2::geom_density(
        ggplot2::aes(fill = .data$assembly),
        alpha = 0.5
      )
    } else {
      p <- p + ggplot2::geom_density(fill = "#56B4E9", alpha = 0.7)
    }
  } else if (plot_type == "histogram") {
    if (by_assembly) {
      p <- p + ggplot2::geom_histogram(
        ggplot2::aes(fill = .data$assembly),
        position = "identity",
        alpha = 0.5,
        bins = 30
      )
    } else {
      p <- p + ggplot2::geom_histogram(fill = "#56B4E9", alpha = 0.7, bins = 30)
    }
  }

  # Add threshold line if specified
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_vline(
      xintercept = threshold,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    )
  }

  # Facet if by_assembly
  if (by_assembly) {
    p <- p + ggplot2::facet_wrap(~ assembly)
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(
      x = score_column,
      y = if (plot_type == "density") "Density" else "Count",
      title = paste("Distribution of", score_column)
    ) +
    ggplot2::theme_minimal()

  p
}

#' Extract scores from protein collection
#'
#' @param proteins A protein_collection object
#' @param score_column Name of the score column
#' @return Data frame with columns: assembly, score
#' @keywords internal
extract_scores <- function(proteins, score_column) {
  results <- lapply(proteins$assemblies, function(ps) {
    if (score_column %in% names(ps$proteins)) {
      data.frame(
        assembly = ps$assembly_name,
        score = ps$proteins[[score_column]],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })

  # Remove NULLs and combine
  results <- results[!vapply(results, is.null, logical(1))]

  if (length(results) == 0) {
    return(NULL)
  }

  do.call(rbind, results)
}
