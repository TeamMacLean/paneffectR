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
#' @section File naming conventions:
#' Each FASTA file in `fasta_dir` is treated as one assembly. The assembly name
#' is the file name without its extension (and without a trailing `.gz`), so
#' `isolate_A.faa` and `isolate_A.faa.gz` both give the assembly `isolate_A`.
#'
#' Optional companion files are matched by assembly name:
#'
#' * `<assembly>_scored.csv` in `score_dir` - effector prediction scores
#'   (see below). A warning is issued for every assembly with no score file.
#' * `<assembly>_name_mapping.csv` in `fasta_dir` - omnieff name mapping
#'   (recorded in `source_files`, not merged).
#'
#' @section FASTA requirements:
#' Files must contain protein (amino acid) sequences. Nucleotide input is
#' rejected with an error. Headers are split at the first whitespace: the first
#' token becomes `protein_id` and the remainder (if any) is kept in a
#' `description` column. For example `>XP_001 cytochrome P450 [P. oryzae]`
#' gives `protein_id = "XP_001"` and `description = "cytochrome P450 [P. oryzae]"`.
#' Protein IDs must be unique across all assemblies. Trailing stop codons (`*`)
#' are removed and sequences are upper-cased.
#'
#' @section Score file format:
#' A score file is a CSV with one row per protein and at least these columns:
#'
#' * `protein_id` - must match the FASTA IDs (first header token) exactly
#' * `custom_score` - numeric effector score used by `score_threshold` in
#'   [build_pa_matrix()] and by [plot_scores()]
#'
#' The omnieff pipeline also writes `assembly`, `score_rank`, and many
#' per-tool columns; all extra columns are kept. The value `RESULTS_MISSING`
#' (and empty strings) are read as `NA`. A warning is issued if none of the
#' IDs in a score file match the FASTA.
#'
#' @param fasta_dir Character. Directory containing FASTA files.
#' @param score_dir Character. Optional directory containing score CSV files.
#' @param pattern Character. Glob pattern for FASTA files (default: "*.faa").
#'   Use e.g. `"*.faa.gz"` for gzipped files or `"*.fa*"` for mixed extensions.
#' @param name_mapping Logical. Whether to record name mapping files.
#' @param quiet Logical. If `FALSE` (default), print a one-line summary of what
#'   was loaded. Warnings are shown regardless.
#'
#' @return A `protein_collection` object.
#' @export
#'
#' @examples
#' # Bundled example data: three assemblies with scores
#' data_dir <- system.file("testdata", package = "paneffectR")
#' proteins <- load_proteins(fasta_dir = data_dir, score_dir = data_dir)
#' proteins
#'
#' \dontrun{
#' # Load from omnieff output
#' proteins <- load_proteins(
#'   fasta_dir = "reformatted/",
#'   score_dir = "scored/"
#' )
#'
#' # Load gzipped FASTAs without scores
#' proteins <- load_proteins(fasta_dir = "my_fastas/", pattern = "*.faa.gz")
#' }
load_proteins <- function(fasta_dir,
                          score_dir = NULL,
                          pattern = "*.faa",
                          name_mapping = TRUE,
                          quiet = FALSE) {
  # Validate fasta_dir
  if (!is.character(fasta_dir) || length(fasta_dir) != 1) {
    cli::cli_abort("{.arg fasta_dir} must be a single character string")
  }
  fasta_dir <- path.expand(fasta_dir)

  if (!dir.exists(fasta_dir)) {
    cli::cli_abort("Directory not found: {.path {fasta_dir}}")
  }

  if (!is.null(score_dir)) {
    if (!is.character(score_dir) || length(score_dir) != 1) {
      cli::cli_abort("{.arg score_dir} must be a single character string or NULL")
    }
    score_dir <- path.expand(score_dir)
    if (!dir.exists(score_dir)) {
      cli::cli_abort("Score directory not found: {.path {score_dir}}")
    }
  }

  # Find FASTA files matching pattern
  fasta_files <- Sys.glob(file.path(fasta_dir, pattern))

  if (length(fasta_files) == 0) {
    cli::cli_abort("No files matching pattern {.val {pattern}} found in {.path {fasta_dir}}")
  }

  # Bookkeeping for summary / warnings
  missing_scores <- character()
  n_scored <- 0L
  n_score_matched <- 0L
  n_with_description <- 0L
  empty_assemblies <- character()

  # Process each FASTA file
  protein_sets <- lapply(fasta_files, function(fasta_path) {
    # Extract assembly name from filename (stem without extension, ignoring .gz)
    assembly_name <- fasta_assembly_name(fasta_path)

    # Load FASTA
    proteins <- load_fasta(fasta_path)

    if (nrow(proteins) == 0) {
      empty_assemblies <<- c(empty_assemblies, assembly_name)
    }
    n_with_description <<- n_with_description + sum(!is.na(proteins$description))

    # Track source files
    source_files <- list(fasta = fasta_path)

    # Try to load scores if score_dir provided
    if (!is.null(score_dir)) {
      score_file <- file.path(score_dir, paste0(assembly_name, "_scored.csv"))
      if (file.exists(score_file)) {
        scores <- load_scores(score_file)
        source_files$scores <- score_file

        # Avoid duplicate columns on join (e.g. an 'assembly' or 'sequence' column)
        clash <- setdiff(intersect(names(scores), names(proteins)), "protein_id")
        if (length(clash) > 0) {
          scores <- scores[, setdiff(names(scores), clash), drop = FALSE]
        }

        n_matched <- sum(proteins$protein_id %in% scores$protein_id)
        if (nrow(proteins) > 0 && n_matched == 0) {
          cli::cli_warn(c(
            "No protein IDs in {.path {score_file}} match the FASTA IDs of assembly {.val {assembly_name}}",
            "i" = "All scores for this assembly will be {.val NA}.",
            "i" = "FASTA IDs look like {.val {utils::head(proteins$protein_id, 2)}}; score IDs look like {.val {utils::head(scores$protein_id, 2)}}."
          ))
        }
        n_scored <<- n_scored + 1L
        n_score_matched <<- n_score_matched + n_matched

        # Join scores to proteins by protein_id
        proteins <- dplyr::left_join(proteins, scores, by = "protein_id")
      } else {
        missing_scores <<- c(missing_scores, assembly_name)
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
  collection <- new_protein_collection(protein_sets)

  # Warnings about surprising inputs
  if (length(empty_assemblies) > 0) {
    cli::cli_warn(
      "{length(empty_assemblies)} FASTA file{?s} contain{?s/} no sequences: {.val {empty_assemblies}}"
    )
  }

  if (length(missing_scores) > 0) {
    cli::cli_warn(c(
      "No score file found for {length(missing_scores)} of {length(fasta_files)} assembl{?y/ies}: {.val {missing_scores}}",
      "i" = "Expected {.file <assembly>_scored.csv} in {.path {score_dir}}. Scores for these assemblies will be {.val NA}."
    ))
  }

  if (length(fasta_files) == 1) {
    cli::cli_warn(c(
      "Only one assembly was loaded ({.val {collection$summary$assembly_name}})",
      "i" = "Comparative analysis needs at least two assemblies; clustering a single assembly yields only singletons."
    ))
  }

  # One-line summary
  if (!quiet) {
    total <- sum(collection$summary$n_proteins)
    msg <- c(
      "v" = "Loaded {length(fasta_files)} assembl{?y/ies} ({total} protein{?s}) from {.path {fasta_dir}}"
    )
    if (!is.null(score_dir)) {
      msg <- c(msg, "i" = "Scores: found for {n_scored} of {length(fasta_files)} assembl{?y/ies}, matched {n_score_matched} of {total} protein{?s}")
    }
    if (n_with_description > 0) {
      msg <- c(msg, "i" = "FASTA headers trimmed to first token for {n_with_description} protein{?s}; the rest of each header is kept in the {.field description} column")
    }
    cli::cli_inform(msg)
  }

  collection
}

#' Derive an assembly name from a FASTA file path
#'
#' Strips the directory, the extension, and any `.gz`/`.bz2`/`.xz` suffix.
#'
#' @param path Character path to a FASTA file.
#' @return Character assembly name.
#' @keywords internal
fasta_assembly_name <- function(path) {
  tools::file_path_sans_ext(basename(path), compression = TRUE)
}

#' Load a single FASTA file
#'
#' Parses a protein FASTA file into a tibble of protein IDs and sequences.
#' Headers are split at the first whitespace: the first token becomes
#' `protein_id` and the remainder is kept as `description` (`NA` if the header
#' had no description). Trailing stop codons (`*`) are removed and sequences
#' are upper-cased. Gzipped files are supported.
#'
#' @param path Character. Path to FASTA file.
#'
#' @return A tibble with columns: protein_id, sequence, description.
#' @export
#'
#' @examples
#' fasta <- system.file("testdata", "assembly1.faa", package = "paneffectR")
#' seqs <- load_fasta(fasta)
#' head(seqs)
load_fasta <- function(path) {
  # Validate input
  if (!is.character(path) || length(path) != 1) {
    cli::cli_abort("{.arg path} must be a single character string")
  }
  path <- path.expand(path)

  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  # Read FASTA using Biostrings
  seqs <- Biostrings::readAAStringSet(path)

  headers <- names(seqs)
  if (is.null(headers)) headers <- character()
  headers <- trimws(headers)

  # Split header at first whitespace: ID is the first token, rest is description
  has_desc <- grepl("\\s", headers)
  protein_id <- sub("\\s.*$", "", headers)
  description <- rep(NA_character_, length(headers))
  description[has_desc] <- trimws(sub("^\\S+\\s+", "", headers[has_desc]))

  # Normalise sequences: upper-case, strip trailing stop codon(s)
  sequence <- toupper(as.character(seqs))
  sequence <- sub("\\*+$", "", sequence)

  # Reject nucleotide input early (DIAMOND hangs on it)
  nuc_frac <- nucleotide_fraction(sequence)
  if (!is.na(nuc_frac) && nuc_frac >= 0.95) {
    cli::cli_abort(c(
      "{.path {path}} looks like nucleotide sequence, not protein",
      "x" = "{round(nuc_frac * 100)}% of residues are A/C/G/T/U/N.",
      "i" = "paneffectR needs amino acid FASTA files. Translate the sequences first."
    ))
  }

  tibble::tibble(
    protein_id = unname(protein_id),
    sequence = unname(sequence),
    description = description
  )
}

#' Fraction of residues that are nucleotide letters
#'
#' Used to detect nucleotide FASTA files mislabelled as protein. Looks at up to
#' `max_seqs` sequences.
#'
#' @param sequences Character vector of sequences.
#' @param max_seqs Maximum number of sequences to inspect.
#' @return Numeric fraction in `[0, 1]`, or `NA` if there are no residues.
#' @keywords internal
nucleotide_fraction <- function(sequences, max_seqs = 1000L) {
  if (length(sequences) == 0) return(NA_real_)
  sample <- paste(toupper(utils::head(sequences, max_seqs)), collapse = "")
  total <- nchar(sample)
  if (total == 0) return(NA_real_)
  nuc <- nchar(gsub("[^ACGTUN]", "", sample))
  nuc / total
}

#' Load effector scores from CSV
#'
#' Parses a scored CSV file, such as those produced by the omnieff pipeline.
#' The file must have `protein_id` and `custom_score` columns; see the
#' "Score file format" section of [load_proteins()] for details.
#'
#' @param path Character. Path to scored CSV file.
#'
#' @return A tibble with score columns.
#' @export
#'
#' @examples
#' csv <- system.file("testdata", "assembly1_scored.csv", package = "paneffectR")
#' scores <- load_scores(csv)
#' scores[, c("protein_id", "custom_score", "score_rank")]
load_scores <- function(path) {
  # Validate input
  if (!is.character(path) || length(path) != 1) {
    cli::cli_abort("{.arg path} must be a single character string")
  }
  path <- path.expand(path)

  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  # Peek at the header to validate columns before choosing parsers
  header <- names(readr::read_csv(
    path,
    n_max = 0,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ))

  required_cols <- c("protein_id", "custom_score")
  missing_cols <- setdiff(required_cols, header)
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "Score file {.path {path}} is missing required column{?s}: {.val {missing_cols}}",
      "i" = "Score files need at least {.field protein_id} and {.field custom_score}; see {.help paneffectR::load_proteins}."
    ))
  }

  # Be explicit about key columns that are present; guess the rest
  known <- list(
    assembly = readr::col_character(),
    protein_id = readr::col_character(),
    score_rank = readr::col_integer(),
    custom_score = readr::col_double()
  )
  col_spec <- do.call(
    readr::cols,
    c(known[intersect(names(known), header)], list(.default = readr::col_guess()))
  )

  # RESULTS_MISSING values will become NA automatically
  scores <- readr::read_csv(
    path,
    col_types = col_spec,
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
#'
#' @examples
#' csv <- system.file("testdata", "assembly1_name_mapping.csv", package = "paneffectR")
#' mapping <- load_name_mapping(csv)
#' head(mapping)
load_name_mapping <- function(path) {
  # Validate input
  if (!is.character(path) || length(path) != 1) {
    cli::cli_abort("{.arg path} must be a single character string")
  }
  path <- path.expand(path)

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
