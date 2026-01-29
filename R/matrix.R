#' Presence/Absence Matrix Functions
#'
#' @name matrix-functions
#' @keywords internal
NULL

#' Build presence/absence matrix from clustering results
#'
#' Creates a matrix with orthogroups as rows and assemblies as columns.
#' Can produce binary (0/1), score-based, or count matrices.
#'
#' @param orthogroups An `orthogroup_result` object.
#' @param proteins Optional `protein_collection` with score data.
#' @param type Character. Matrix type: "binary", "count", or "score".
#' @param score_column Character. Column name for scores (default: "custom_score").
#' @param score_threshold Numeric. Only include proteins with score >= threshold.
#' @param score_aggregation Character. How to aggregate scores: "max", "mean", or "sum".
#'
#' @return A `pa_matrix` object.
#' @export
#'
#' @examples
#' \dontrun{
#' pa <- build_pa_matrix(clusters)
#' pa <- build_pa_matrix(clusters, score_threshold = 5.0)
#' }
build_pa_matrix <- function(orthogroups,
                            proteins = NULL,
                            type = "binary",
                            score_column = "custom_score",
                            score_threshold = NULL,
                            score_aggregation = "max") {
  # Validate input is orthogroup_result
  if (!inherits(orthogroups, "orthogroup_result")) {
    cli::cli_abort("{.arg orthogroups} must be an {.cls orthogroup_result} object")
  }

  # Validate type
  type <- match.arg(type, c("binary", "count", "score"))

  # Score type not yet implemented
  if (type == "score") {
    cli::cli_abort("Type {.val {type}} is not yet implemented")
  }

  # Extract orthogroups tibble
  og_tibble <- orthogroups$orthogroups

  # Get unique orthogroup IDs and assembly names
  og_ids <- unique(og_tibble$orthogroup_id)
  assemblies <- unique(og_tibble$assembly)

  # Sort for consistent ordering
  og_ids <- sort(og_ids)
  assemblies <- sort(assemblies)

  # Create empty matrix initialized to 0
  mat <- matrix(
    0L,
    nrow = length(og_ids),
    ncol = length(assemblies),
    dimnames = list(og_ids, assemblies)
  )

  # Fill in matrix based on type
  if (type == "binary") {
    # Binary: 1 if any protein from assembly in orthogroup
    for (i in seq_len(nrow(og_tibble))) {
      og <- og_tibble$orthogroup_id[i]
      asm <- og_tibble$assembly[i]
      mat[og, asm] <- 1L
    }
  } else if (type == "count") {
    # Count: increment for each protein (detects paralogs)
    for (i in seq_len(nrow(og_tibble))) {
      og <- og_tibble$orthogroup_id[i]
      asm <- og_tibble$assembly[i]
      mat[og, asm] <- mat[og, asm] + 1L
    }
  }

  # Build orthogroups metadata tibble (orthogroup_id, size)
  og_sizes <- og_tibble |>
    dplyr::count(.data$orthogroup_id, name = "size") |>
    dplyr::arrange(.data$orthogroup_id)

  # Build assemblies metadata tibble (assembly_name, n_orthogroups)
  # Count how many orthogroups each assembly is present in
  asm_counts <- og_tibble |>
    dplyr::distinct(.data$orthogroup_id, .data$assembly) |>
    dplyr::count(.data$assembly, name = "n_orthogroups") |>
    dplyr::rename(assembly_name = "assembly") |>
    dplyr::arrange(.data$assembly_name)

  # Create pa_matrix object
  new_pa_matrix(
    matrix = mat,
    orthogroups = og_sizes,
    assemblies = asm_counts,
    type = type,
    threshold = score_threshold
  )
}

#' Filter presence/absence matrix by score
#'
#' Removes orthogroups where no member meets the score threshold.
#'
#' @param pa A `pa_matrix` object.
#' @param proteins A `protein_collection` with score data.
#' @param threshold Numeric. Minimum score to retain an orthogroup.
#' @param score_column Character. Column name for scores.
#'
#' @return A filtered `pa_matrix` object.
#' @export
filter_by_score <- function(pa,
                            proteins,
                            threshold,
                            score_column = "custom_score") {
  cli::cli_abort("Not implemented")
}
