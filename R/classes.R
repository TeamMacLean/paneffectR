#' S3 Class Constructors for paneffectR
#'
#' @name paneffectR-classes
#' @keywords internal
NULL

# protein_set class --------------------------------------------------------

#' Create a protein_set object
#'
#' Represents proteins from a single assembly.
#'
#' @param assembly_name Character. Unique identifier for the assembly.
#' @param proteins A tibble with columns: protein_id, sequence, and optionally
#'   score, rank, and other metadata.
#' @param metadata List. Optional assembly-level metadata.
#' @param source_files List. Paths to original .faa and .csv files.
#'
#' @return A `protein_set` S3 object.
#' @export
new_protein_set <- function(assembly_name,
                            proteins,
                            metadata = list(),
                            source_files = list()) {
  # Validate assembly_name
  if (!is.character(assembly_name) || length(assembly_name) != 1 || nchar(assembly_name) == 0) {
    cli::cli_abort("{.arg assembly_name} must be a single non-empty character string")
  }

 # Validate proteins is a tibble with required columns
  if (!tibble::is_tibble(proteins)) {
    cli::cli_abort("{.arg proteins} must be a tibble")
  }

  required_cols <- c("protein_id", "sequence")
  missing_cols <- setdiff(required_cols, names(proteins))
  if (length(missing_cols) > 0) {
    cli::cli_abort("{.arg proteins} must have columns: {.val {required_cols}}")
  }

  # Validate protein_id uniqueness
  if (anyDuplicated(proteins$protein_id)) {
    cli::cli_abort("Duplicate protein IDs found in {.arg proteins}")
  }

  # Validate metadata and source_files are lists
  if (!is.list(metadata)) {
    cli::cli_abort("{.arg metadata} must be a list")
  }
  if (!is.list(source_files)) {
    cli::cli_abort("{.arg source_files} must be a list")
  }

  # Construct object
  structure(
    list(
      assembly_name = assembly_name,
      proteins = proteins,
      metadata = metadata,
      source_files = source_files
    ),
    class = "protein_set"
  )
}

#' @export
print.protein_set <- function(x, ...) {
  n_proteins <- nrow(x$proteins)
  has_scores <- "custom_score" %in% names(x$proteins)

  cat("-- protein_set:", x$assembly_name, "--\n")
  cat(n_proteins, " protein", if (n_proteins != 1) "s", "\n", sep = "")

  if (has_scores) {
    cat("Scores: present\n")
  }

  invisible(x)
}

# protein_collection class -------------------------------------------------

#' Create a protein_collection object
#'
#' Container for multiple assemblies.
#'
#' @param assemblies List of protein_set objects.
#'
#' @return A `protein_collection` S3 object.
#' @export
new_protein_collection <- function(assemblies) {
  # Validate input is a non-empty list
  if (!is.list(assemblies) || length(assemblies) == 0) {
    cli::cli_abort("{.arg assemblies} must be a non-empty list")
  }

  # Validate all items are protein_set objects
  for (i in seq_along(assemblies)) {
    if (!inherits(assemblies[[i]], "protein_set")) {
      cli::cli_abort("Item {i} in {.arg assemblies} is not a {.cls protein_set}")
    }
  }

  # Extract assembly names
  assembly_names <- vapply(assemblies, function(x) x$assembly_name, character(1))

  # Check for duplicate assembly names
  if (anyDuplicated(assembly_names)) {
    dups <- assembly_names[duplicated(assembly_names)]
    cli::cli_abort("Duplicate assembly names: {.val {unique(dups)}}")
  }

  # Check for globally unique protein_ids
  all_protein_ids <- unlist(lapply(assemblies, function(x) x$proteins$protein_id))
  if (anyDuplicated(all_protein_ids)) {
    dups <- all_protein_ids[duplicated(all_protein_ids)]
    cli::cli_abort("Duplicate protein IDs across assemblies: {.val {head(unique(dups), 5)}}")
  }

  # Name the list by assembly names
  names(assemblies) <- assembly_names

  # Compute summary
  summary_df <- tibble::tibble(
    assembly_name = assembly_names,
    n_proteins = vapply(assemblies, function(x) nrow(x$proteins), integer(1)),
    has_scores = vapply(assemblies, function(x) "custom_score" %in% names(x$proteins), logical(1))
  )

  # Construct object
  structure(
    list(
      assemblies = assemblies,
      n_assemblies = length(assemblies),
      summary = summary_df
    ),
    class = "protein_collection"
  )
}

#' @export
print.protein_collection <- function(x, ...) {
  total_proteins <- sum(x$summary$n_proteins)
  n_asm <- x$n_assemblies

  cat("-- protein_collection --\n")
  cat(n_asm, " assembl", if (n_asm == 1) "y" else "ies",
      ", ", total_proteins, " total protein",
      if (total_proteins != 1) "s", "\n", sep = "")
  cat("\n")

  # Print summary table
  print(x$summary, n = min(10, nrow(x$summary)))

  invisible(x)
}

# orthogroup_result class --------------------------------------------------

#' Create an orthogroup_result object
#'
#' Result of clustering proteins across assemblies.
#'
#' @param orthogroups A tibble with columns: orthogroup_id, assembly, protein_id.
#' @param method Character. Clustering method used.
#' @param parameters List. Clustering parameters used.
#' @param singletons A tibble. Proteins not assigned to any group.
#' @param stats A tibble. Orthogroup sizes and assembly coverage.
#'
#' @return An `orthogroup_result` S3 object.
#' @export
new_orthogroup_result <- function(orthogroups,
                                  method,
                                  parameters = list(),
                                  singletons = NULL,
                                  stats = NULL) {
  cli::cli_abort("orthogroup_result is not yet implemented (Phase 2)")
}

#' @export
print.orthogroup_result <- function(x, ...) {
  cli::cli_abort("orthogroup_result is not yet implemented (Phase 2)")
}

# pa_matrix class ----------------------------------------------------------

#' Create a pa_matrix object
#'
#' Presence/absence matrix with metadata.
#'
#' @param matrix Matrix. Rows = orthogroups, cols = assemblies (0/1 or scores).
#' @param orthogroups A tibble. Orthogroup metadata.
#' @param assemblies A tibble. Assembly metadata.
#' @param type Character. One of "binary", "score", or "count".
#' @param threshold Numeric. Score threshold used (if applicable).
#'
#' @return A `pa_matrix` S3 object.
#' @export
new_pa_matrix <- function(matrix,
                          orthogroups,
                          assemblies,
                          type = "binary",
                          threshold = NULL) {
  cli::cli_abort("pa_matrix is not yet implemented (Phase 3)")
}

#' @export
print.pa_matrix <- function(x, ...) {
  cli::cli_abort("pa_matrix is not yet implemented (Phase 3)")
}
