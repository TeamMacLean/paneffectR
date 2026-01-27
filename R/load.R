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
  cli::cli_abort("Not implemented")
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
  cli::cli_abort("Not implemented")
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
  cli::cli_abort("Not implemented")
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
  cli::cli_abort("Not implemented")
}
