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
  # Validate orthogroups is a tibble

  if (!tibble::is_tibble(orthogroups)) {
    cli::cli_abort("{.arg orthogroups} must be a tibble")
  }

  # Validate required columns in orthogroups
 required_cols <- c("orthogroup_id", "assembly", "protein_id")
  missing_cols <- setdiff(required_cols, names(orthogroups))
  if (length(missing_cols) > 0) {
    cli::cli_abort("{.arg orthogroups} must have columns: {.val {missing_cols}}")
  }

  # Validate method is single character string
  if (!is.character(method) || length(method) != 1) {
    cli::cli_abort("{.arg method} must be a single character string")
  }

  # Validate parameters is list
  if (!is.list(parameters)) {
    cli::cli_abort("{.arg parameters} must be a list")
  }

  # Validate singletons if provided
  if (!is.null(singletons)) {
    if (!tibble::is_tibble(singletons)) {
      cli::cli_abort("{.arg singletons} must be a tibble")
    }
    singleton_required <- c("assembly", "protein_id")
    singleton_missing <- setdiff(singleton_required, names(singletons))
    if (length(singleton_missing) > 0) {
      cli::cli_abort("{.arg singletons} must have columns: {.val {singleton_required}}")
    }
  }

  # Auto-compute stats if not provided
  if (is.null(stats)) {
    n_orthogroups <- length(unique(orthogroups$orthogroup_id))
    n_singletons <- if (is.null(singletons)) 0L else nrow(singletons)
    n_proteins_clustered <- nrow(orthogroups)
    n_assemblies <- length(unique(orthogroups$assembly))

    stats <- tibble::tibble(
      n_orthogroups = n_orthogroups,
      n_singletons = n_singletons,
      n_proteins_clustered = n_proteins_clustered,
      n_assemblies = n_assemblies
    )
  }

  # Construct object
  structure(
    list(
      orthogroups = orthogroups,
      method = method,
      parameters = parameters,
      singletons = singletons,
      stats = stats
    ),
    class = "orthogroup_result"
  )
}

#' @export
print.orthogroup_result <- function(x, ...) {
  n_og <- x$stats$n_orthogroups
  n_sing <- x$stats$n_singletons

  cat("-- orthogroup_result (", x$method, ") --\n", sep = "")
  cat(n_og, " orthogroup", if (n_og != 1) "s", "\n", sep = "")

  if (n_sing > 0) {
    cat(n_sing, " singleton", if (n_sing != 1) "s", "\n", sep = "")
  }

  invisible(x)
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
  # Validate matrix is a matrix
  if (!is.matrix(matrix)) {
    cli::cli_abort("{.arg matrix} must be a matrix")
  }

  # Validate orthogroups is a tibble with required columns
  if (!tibble::is_tibble(orthogroups)) {
    cli::cli_abort("{.arg orthogroups} must be a tibble")
  }
  if (!"orthogroup_id" %in% names(orthogroups)) {
    cli::cli_abort("{.arg orthogroups} must have column: {.val orthogroup_id}")
  }

  # Validate assemblies is a tibble with required columns
  if (!tibble::is_tibble(assemblies)) {
    cli::cli_abort("{.arg assemblies} must be a tibble")
  }
  if (!"assembly_name" %in% names(assemblies)) {
    cli::cli_abort("{.arg assemblies} must have column: {.val assembly_name}")
  }

  # Validate matrix dimensions match metadata
  if (nrow(matrix) != nrow(orthogroups)) {
    cli::cli_abort(
      "Matrix row count ({nrow(matrix)}) does not match orthogroups count ({nrow(orthogroups)})"
    )
  }
  if (ncol(matrix) != nrow(assemblies)) {
    cli::cli_abort(
      "Matrix column count ({ncol(matrix)}) does not match assemblies count ({nrow(assemblies)})"
    )
  }

  # Validate type is one of allowed values
  allowed_types <- c("binary", "count", "score")
  if (!is.character(type) || length(type) != 1 || !type %in% allowed_types) {
    cli::cli_abort("{.arg type} must be one of: {.val {allowed_types}}")
  }

  # Validate threshold is NULL or numeric
  if (!is.null(threshold) && !is.numeric(threshold)) {
    cli::cli_abort("{.arg threshold} must be NULL or numeric")
  }

  # Construct object
  structure(
    list(
      matrix = matrix,
      orthogroups = orthogroups,
      assemblies = assemblies,
      type = type,
      threshold = threshold
    ),
    class = "pa_matrix"
  )
}

#' @export
print.pa_matrix <- function(x, ...) {
  n_og <- nrow(x$matrix)
  n_asm <- ncol(x$matrix)
  total_cells <- n_og * n_asm

  # Calculate sparsity (percentage of zeros)
  n_zeros <- sum(x$matrix == 0)
  sparsity_pct <- round(100 * n_zeros / total_cells, 1)

  cat("-- pa_matrix (", x$type, ") --\n", sep = "")
  cat(n_og, " orthogroup", if (n_og != 1) "s", " x ",
      n_asm, " assembl", if (n_asm == 1) "y" else "ies", "\n", sep = "")
  cat("Sparsity: ", sparsity_pct, "%\n", sep = "")

  if (!is.null(x$threshold)) {
    cat("Threshold: ", x$threshold, "\n", sep = "")
  }

  invisible(x)
}

#' Subset a pa_matrix object
#'
#' @param x A `pa_matrix` object.
#' @param i Row indices (orthogroups).
#' @param j Column indices (assemblies).
#' @param drop Ignored (always FALSE to preserve pa_matrix structure).
#' @param ... Additional arguments (ignored).
#'
#' @return A subsetted `pa_matrix` object.
#' @export
`[.pa_matrix` <- function(x, i, j, drop = FALSE, ...) {
  # Handle missing indices
  if (missing(i)) i <- seq_len(nrow(x$matrix))
  if (missing(j)) j <- seq_len(ncol(x$matrix))

  # Subset the matrix (drop = FALSE to keep matrix structure)
  new_mat <- x$matrix[i, j, drop = FALSE]

  # Get the row and column names of the subsetted matrix
  new_og_ids <- rownames(new_mat)
  new_asm_names <- colnames(new_mat)

  # Subset orthogroups metadata to match
  new_orthogroups <- x$orthogroups |>
    dplyr::filter(.data$orthogroup_id %in% new_og_ids) |>
    dplyr::arrange(match(.data$orthogroup_id, new_og_ids))

  # Subset assemblies metadata to match
  new_assemblies <- x$assemblies |>
    dplyr::filter(.data$assembly_name %in% new_asm_names) |>
    dplyr::arrange(match(.data$assembly_name, new_asm_names))

  # Return new pa_matrix
  new_pa_matrix(
    matrix = new_mat,
    orthogroups = new_orthogroups,
    assemblies = new_assemblies,
    type = x$type,
    threshold = x$threshold
  )
}
#' Convert pa_matrix to data.frame
#'
#' @param x A `pa_matrix` object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param format Character. "long" (default) or "wide".
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame in long or wide format.
#' @export
as.data.frame.pa_matrix <- function(x, row.names = NULL, optional = FALSE,
                                    format = "long", ...) {
  format <- match.arg(format, c("long", "wide"))

  if (format == "wide") {
    # Wide format: orthogroup_id + assembly columns
    df <- as.data.frame(x$matrix)
    df$orthogroup_id <- rownames(x$matrix)
    # Reorder columns to put orthogroup_id first
    df <- df[, c("orthogroup_id", colnames(x$matrix)), drop = FALSE]
    rownames(df) <- NULL
  } else {
    # Long format: orthogroup_id, assembly, value
    df <- expand.grid(
      orthogroup_id = rownames(x$matrix),
      assembly = colnames(x$matrix),
      stringsAsFactors = FALSE
    )
    df$value <- as.vector(x$matrix)
  }

  df
}
