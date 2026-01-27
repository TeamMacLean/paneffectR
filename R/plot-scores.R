#' Score Distribution Visualization
#'
#' @name plot-scores
#' @keywords internal
NULL

#' Plot effector score distributions
#'
#' Creates density plots or histograms of effector scores per assembly.
#'
#' @param proteins A `protein_collection` object with score data.
#' @param score_column Character. Column name for scores (default: "custom_score").
#' @param type Character. Plot type: "density" or "histogram".
#' @param facet Logical. Whether to facet by assembly (default: FALSE).
#' @param threshold Numeric. Optional threshold line to draw.
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_scores(proteins)
#' plot_scores(proteins, threshold = 5.0)
#' }
plot_scores <- function(proteins,
                        score_column = "custom_score",
                        type = c("density", "histogram"),
                        facet = FALSE,
                        threshold = NULL,
                        ...) {
  cli::cli_abort("Not implemented")
}
