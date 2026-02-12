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

  # Validate score_aggregation
 score_aggregation <- match.arg(score_aggregation, c("max", "mean", "sum"))

  # Score type requires proteins argument
  if (type == "score" && is.null(proteins)) {
    cli::cli_abort("Type {.val score} requires {.arg proteins} argument")
  }

  # score_threshold requires proteins argument
  if (!is.null(score_threshold) && is.null(proteins)) {
    cli::cli_abort("{.arg score_threshold} requires {.arg proteins} argument")
  }

  # Extract orthogroups tibble
  og_tibble <- orthogroups$orthogroups

  # Apply score_threshold filtering if specified
  if (!is.null(score_threshold)) {
    score_lookup <- build_score_lookup(proteins, score_column)

    # Calculate max score per orthogroup
    og_max_scores <- og_tibble |>
      dplyr::mutate(score = score_lookup[.data$protein_id]) |>
      dplyr::group_by(.data$orthogroup_id) |>
      dplyr::summarise(max_score = max(.data$score, na.rm = TRUE), .groups = "drop")

    # Keep only orthogroups meeting threshold
    passing_ogs <- og_max_scores$orthogroup_id[og_max_scores$max_score >= score_threshold]

    # Filter orthogroups tibble
    og_tibble <- og_tibble |>
      dplyr::filter(.data$orthogroup_id %in% passing_ogs)
  }

  # Get unique orthogroup IDs and assembly names
  og_ids <- unique(og_tibble$orthogroup_id)
  asm_names <- unique(og_tibble$assembly)

  # Sort for consistent ordering
  og_ids <- sort(og_ids)
  asm_names <- sort(asm_names)

  # Build matrix based on type
  if (type == "score") {
    # Build protein_id -> score lookup from protein_collection
    score_lookup <- build_score_lookup(proteins, score_column)

    # Create empty matrix initialized to NA (absent = NA for scores)
    mat <- matrix(
      NA_real_,
      nrow = length(og_ids),
      ncol = length(asm_names),
      dimnames = list(og_ids, asm_names)
    )

    # Group by orthogroup and assembly, then aggregate scores
    for (og in og_ids) {
      for (asm in asm_names) {
        # Get protein_ids for this orthogroup + assembly
        protein_ids <- og_tibble$protein_id[
          og_tibble$orthogroup_id == og & og_tibble$assembly == asm
        ]

        if (length(protein_ids) > 0) {
          # Look up scores
          scores <- score_lookup[protein_ids]
          scores <- scores[!is.na(scores)]

          if (length(scores) > 0) {
            # Aggregate scores
            mat[og, asm] <- switch(
              score_aggregation,
              "max" = max(scores),
              "mean" = mean(scores),
              "sum" = sum(scores)
            )
          }
        }
      }
    }
  } else {
    # Binary or count: initialize to 0
    mat <- matrix(
      0L,
      nrow = length(og_ids),
      ncol = length(asm_names),
      dimnames = list(og_ids, asm_names)
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

#' Build protein_id to score lookup from protein_collection
#'
#' @param proteins A `protein_collection` object.
#' @param score_column Character. Column name for scores.
#'
#' @return Named numeric vector (protein_id -> score).
#' @keywords internal
build_score_lookup <- function(proteins, score_column) {
  # Validate proteins is a protein_collection
 if (!inherits(proteins, "protein_collection")) {
    cli::cli_abort("{.arg proteins} must be a {.cls protein_collection} object")
  }

  # Combine all proteins from all assemblies
  all_proteins <- do.call(rbind, lapply(proteins$assemblies, function(ps) {
    ps$proteins
  }))

  # Check score column exists
  if (!score_column %in% names(all_proteins)) {
    cli::cli_abort(
      "Score column {.val {score_column}} not found in proteins. Available columns: {.val {names(all_proteins)}}"
    )
  }

  # Build named vector: protein_id -> score
  scores <- all_proteins[[score_column]]
  names(scores) <- all_proteins$protein_id

  scores
}
