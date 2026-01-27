#' Data Loading Functions
#'
#' @name load-functions
#' @keywords internal
NULL

#' Load proteins from multiple assemblies
#'
#' Main entry point for loading protein data. Discovers and loads FASTA files,
#' optionally with associated score and name mapping files.
#'
#' @param fasta_dir Character. Directory containing FASTA files.
#' @param score_dir Character. Optional directory containing score CSV files.
#' @param pattern Character. Glob pattern for FASTA files (default: "*.faa").
#' @param name_mapping Logical. Whether to load name mapping files.
#'
#' @return A `protein_collection` object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load from omnieff output
#' proteins <- load_proteins(
#'   fasta_dir = "reformatted/",
#'   score_dir = "scored/"
#' )
#'
#' # Load just FASTAs
#' proteins <- load_proteins(fasta_dir = "my_fastas/")
#' }
load_proteins <- function(fasta_dir,
                          score_dir = NULL,
                          pattern = "*.faa",
                          name_mapping = TRUE) {
  # Validate fasta_dir
  if (!is.character(fasta_dir) || length(fasta_dir) != 1) {
    cli::cli_abort("{.arg fasta_dir} must be a single character string")
  }

  if (!dir.exists(fasta_dir)) {
    cli::cli_abort("Directory not found: {.path {fasta_dir}}")
  }

  # Find FASTA files matching pattern
  fasta_files <- Sys.glob(file.path(fasta_dir, pattern))

  if (length(fasta_files) == 0) {
    cli::cli_abort("No files matching pattern {.val {pattern}} found in {.path {fasta_dir}}")
  }

  # Process each FASTA file
  protein_sets <- lapply(fasta_files, function(fasta_path) {
    # Extract assembly name from filename (stem without extension)
    assembly_name <- tools::file_path_sans_ext(basename(fasta_path))

    # Load FASTA
    proteins <- load_fasta(fasta_path)

    # Track source files
    source_files <- list(fasta = fasta_path)

    # Try to load scores if score_dir provided
    if (!is.null(score_dir)) {
      score_file <- file.path(score_dir, paste0(assembly_name, "_scored.csv"))
      if (file.exists(score_file)) {
        scores <- load_scores(score_file)
        source_files$scores <- score_file

        # Join scores to proteins by protein_id
        proteins <- dplyr::left_join(proteins, scores, by = "protein_id")
      }
    }

    # Try to load name mapping if requested
    if (name_mapping) {
      # Look for name mapping in same directory as FASTA
      mapping_file <- file.path(dirname(fasta_path), paste0(assembly_name, "_name_mapping.csv"))
      if (file.exists(mapping_file)) {
        source_files$name_mapping <- mapping_file
        # Name mapping is stored but not merged into proteins tibble
        # (it maps old_name -> new_name, proteins already have new_name)
      }
    }

    # Create protein_set
    new_protein_set(
      assembly_name = assembly_name,
      proteins = proteins,
      source_files = source_files
    )
  })

  # Create protein_collection
  new_protein_collection(protein_sets)
}

#' Load a single FASTA file
#'
#' Parses a FASTA file into a tibble of protein IDs and sequences.
#'
#' @param path Character. Path to FASTA file.
#'
#' @return A tibble with columns: protein_id, sequence.
#' @export
load_fasta <- function(path) {
  # Validate input

  if (!is.character(path) || length(path) != 1) {
    cli::cli_abort("{.arg path} must be a single character string")
  }

  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  # Read FASTA using Biostrings
  seqs <- Biostrings::readAAStringSet(path)

  # Convert to tibble
tibble::tibble(
    protein_id = names(seqs),
    sequence = as.character(seqs)
  )
}

#' Load effector scores from CSV
#'
#' Parses a scored CSV file from the omnieff pipeline.
#'
#' @param path Character. Path to scored CSV file.
#'
#' @return A tibble with score columns.
#' @export
load_scores <- function(path) {
  # Validate input
  if (!is.character(path) || length(path) != 1) {
    cli::cli_abort("{.arg path} must be a single character string")
  }

  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  # Read CSV - let readr infer most types but be explicit about key columns
  # RESULTS_MISSING values will become NA automatically
  scores <- readr::read_csv(
    path,
    col_types = readr::cols(
      assembly = readr::col_character(),
      protein_id = readr::col_character(),
      score_rank = readr::col_integer(),
      custom_score = readr::col_double(),
      .default = readr::col_guess()
    ),
    na = c("", "NA", "RESULTS_MISSING"),
    show_col_types = FALSE
  )

  scores
}

#' Load name mapping from CSV
#'
#' Parses a name mapping CSV that maps original protein IDs to standardized names.
#'
#' @param path Character. Path to name mapping CSV file.
#'
#' @return A tibble with columns: assembly, old_name, new_name, protein_number.
#' @export
load_name_mapping <- function(path) {
  # Validate input
  if (!is.character(path) || length(path) != 1) {
    cli::cli_abort("{.arg path} must be a single character string")
  }

  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  # Read CSV
  mapping <- readr::read_csv(
    path,
    col_types = readr::cols(
      assembly = readr::col_character(),
      old_name = readr::col_character(),
      new_name = readr::col_character(),
      protein_number = readr::col_integer()
    ),
    show_col_types = FALSE
  )

  # Validate required columns
  required_cols <- c("assembly", "old_name", "new_name", "protein_number")
  missing_cols <- setdiff(required_cols, names(mapping))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {.val {missing_cols}}")
  }

  mapping
}
