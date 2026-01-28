#' Clustering functions for paneffectR
#'
#' @name paneffectR-clustering
#' @keywords internal
NULL

# RBH Helper Functions ------------------------------------------------------

#' Compute reciprocal best hits from BLAST results
#'
#' Given a table of BLAST hits with bitscores, identifies pairs where
#' each protein is the other's best hit.
#'
#' @param hits A tibble with columns: qseqid, sseqid, bitscore
#'
#' @return A tibble with columns: protein_a, protein_b (RBH pairs)
#' @keywords internal
compute_rbh <- function(hits) {
  # Handle empty input
  if (nrow(hits) == 0) {
    return(tibble::tibble(protein_a = character(), protein_b = character()))
  }

  # Remove self-hits
  hits <- hits[hits$qseqid != hits$sseqid, ]

  if (nrow(hits) == 0) {
    return(tibble::tibble(protein_a = character(), protein_b = character()))
  }

  # Find best hit for each query (highest bitscore)
  best_hits <- hits |>
    dplyr::group_by(.data$qseqid) |>
    dplyr::slice_max(order_by = .data$bitscore, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Create a lookup: query -> best subject
  best_subject <- stats::setNames(best_hits$sseqid, best_hits$qseqid)

  # Find reciprocal pairs
  # A->B is RBH if B's best hit is A
  rbh_pairs <- best_hits |>
    dplyr::filter(
      .data$qseqid %in% names(best_subject) &
      .data$sseqid %in% names(best_subject) &
      best_subject[.data$sseqid] == .data$qseqid
    )

  if (nrow(rbh_pairs) == 0) {
    return(tibble::tibble(protein_a = character(), protein_b = character()))
  }

  # Deduplicate: keep only one direction (alphabetically smaller first)
  rbh_pairs <- rbh_pairs |>
    dplyr::mutate(
      protein_a = pmin(.data$qseqid, .data$sseqid),
      protein_b = pmax(.data$qseqid, .data$sseqid)
    ) |>
    dplyr::distinct(.data$protein_a, .data$protein_b) |>
    dplyr::select("protein_a", "protein_b")

  rbh_pairs
}

#' Build orthogroups from RBH pairs using connected components
#'
#' Groups proteins into orthogroups based on RBH relationships.
#' Proteins connected through any chain of RBH pairs form one orthogroup.
#'
#' @param rbh_pairs A tibble with columns: protein_a, protein_b
#' @param all_protein_ids Character vector of all protein IDs (for identifying singletons)
#'
#' @return A tibble with columns: orthogroup_id, protein_id
#'   Only proteins in orthogroups are returned (singletons excluded)
#' @keywords internal
build_orthogroups_from_rbh <- function(rbh_pairs, all_protein_ids) {
  # Handle empty RBH pairs
  if (nrow(rbh_pairs) == 0) {
    return(tibble::tibble(orthogroup_id = character(), protein_id = character()))
  }

  # Get all proteins involved in RBH pairs
  proteins_in_rbh <- unique(c(rbh_pairs$protein_a, rbh_pairs$protein_b))

  # Build adjacency list
  adj <- list()
  for (p in proteins_in_rbh) {
    adj[[p]] <- character()
  }

  for (i in seq_len(nrow(rbh_pairs))) {
    a <- rbh_pairs$protein_a[i]
    b <- rbh_pairs$protein_b[i]
    adj[[a]] <- c(adj[[a]], b)
    adj[[b]] <- c(adj[[b]], a)
  }

  # Find connected components using BFS
  visited <- character()
  components <- list()

  for (start in proteins_in_rbh) {
    if (start %in% visited) next

    # BFS from this node
    component <- character()
    queue <- start

    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]

      if (node %in% visited) next
      visited <- c(visited, node)
      component <- c(component, node)

      # Add unvisited neighbors to queue
      neighbors <- adj[[node]]
      queue <- c(queue, setdiff(neighbors, visited))
    }

    components[[length(components) + 1]] <- component
  }

  # Build result tibble with sequential OG IDs
  result_list <- lapply(seq_along(components), function(i) {
    og_id <- sprintf("OG%04d", i)
    tibble::tibble(
      orthogroup_id = og_id,
      protein_id = components[[i]]
    )
  })

  dplyr::bind_rows(result_list)
}

# DIAMOND Backend Functions -------------------------------------------------

#' Write protein_collection to single FASTA file
#'
#' Combines all proteins from all assemblies into a single FASTA file
#' suitable for DIAMOND database creation and searching.
#'
#' @param proteins A protein_collection object
#' @param path Path to write the FASTA file
#'
#' @return Invisibly returns the path
#' @keywords internal
write_combined_fasta <- function(proteins, path) {
  # Collect all proteins from all assemblies
  all_proteins <- lapply(proteins$assemblies, function(ps) {
    tibble::tibble(
      protein_id = ps$proteins$protein_id,
      sequence = ps$proteins$sequence
    )
  })
  combined <- dplyr::bind_rows(all_proteins)

  # Write FASTA format
  lines <- character(nrow(combined) * 2)
  for (i in seq_len(nrow(combined))) {
    lines[(i - 1) * 2 + 1] <- paste0(">", combined$protein_id[i])
    lines[(i - 1) * 2 + 2] <- combined$sequence[i]
  }

  writeLines(lines, path)
  invisible(path)
}

#' Parse DIAMOND tabular output
#'
#' Reads DIAMOND blastp output in tabular format and filters by
#' identity, coverage, and e-value thresholds.
#'
#' @param path Path to DIAMOND output file
#' @param min_identity Minimum percent identity threshold
#' @param min_coverage Minimum query coverage threshold
#' @param evalue Maximum e-value threshold
#'
#' @return A tibble with columns: qseqid, sseqid, pident, length, qlen, evalue, bitscore, qcov
#' @keywords internal
parse_diamond_hits <- function(path, min_identity, min_coverage, evalue) {
  # Define empty result structure
  empty_result <- tibble::tibble(
    qseqid = character(),
    sseqid = character(),
    pident = numeric(),
    length = integer(),
    qlen = integer(),
    evalue = numeric(),
    bitscore = numeric(),
    qcov = numeric()
  )

  # Check if file is empty
  if (file.info(path)$size == 0) {
    return(empty_result)
  }

  # Read tabular output
  # Format: qseqid sseqid pident length qlen evalue bitscore
  hits <- readr::read_tsv(
    path,
    col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore"),
    col_types = readr::cols(
      qseqid = readr::col_character(),
      sseqid = readr::col_character(),
      pident = readr::col_double(),
      length = readr::col_integer(),
      qlen = readr::col_integer(),
      evalue = readr::col_double(),
      bitscore = readr::col_double()
    ),
    show_col_types = FALSE
  )

  if (nrow(hits) == 0) {
    return(empty_result)
  }

  # Compute query coverage
  hits <- hits |>
    dplyr::mutate(qcov = (.data$length / .data$qlen) * 100)

  # Apply filters
  # Use .env$ to explicitly reference function parameters to avoid column name conflicts
  hits <- hits |>
    dplyr::filter(
      .data$pident >= .env$min_identity,
      .data$qcov >= .env$min_coverage,
      .data$evalue <= .env$evalue
    )

  hits
}

#' Run DIAMOND reciprocal best hits clustering
#'
#' Performs all-vs-all DIAMOND blastp search and identifies reciprocal
#' best hits to build orthogroups.
#'
#' @param proteins A protein_collection object
#' @param min_identity Minimum percent identity (default 30)
#' @param min_coverage Minimum query coverage (default 50)
#' @param evalue E-value threshold (default 1e-5)
#' @param threads Number of threads (default: auto-detect)
#' @param mode "fast" or "thorough" (default "fast")
#' @param tool_path Optional explicit path to diamond binary
#' @param keep_temp Keep temporary files for debugging (default FALSE)
#'
#' @return An orthogroup_result object
#' @keywords internal
run_diamond_rbh <- function(proteins,
                            min_identity = 30,
                            min_coverage = 50,
                            evalue = 1e-5,
                            threads = NULL,
                            mode = "fast",
                            tool_path = NULL,
                            keep_temp = FALSE) {
  # Find DIAMOND binary
  diamond_path <- if (!is.null(tool_path)) {
    tool_path
  } else if (file.exists("./this_project_env/bin/diamond")) {
    "./this_project_env/bin/diamond"
  } else {
    Sys.which("diamond")
  }

  if (!nzchar(diamond_path) || !file.exists(diamond_path)) {
    cli::cli_abort("DIAMOND not found. Install with: mamba install -c bioconda diamond")
  }

  # Set threads if not specified
  if (is.null(threads)) {
    threads <- parallel::detectCores()
  }

  # Create temp directory

  tmp_dir <- tempfile(pattern = "diamond_rbh_")
  dir.create(tmp_dir)
  if (!keep_temp) {
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  }

  # File paths
  fasta_path <- file.path(tmp_dir, "combined.faa")
  db_path <- file.path(tmp_dir, "db")
  hits_path <- file.path(tmp_dir, "hits.tsv")

  # Write combined FASTA
  write_combined_fasta(proteins, fasta_path)

  # Build database
  makedb_args <- c(
    "makedb",
    "--in", fasta_path,
    "-d", db_path,
    "--quiet"
  )

  makedb_result <- processx::run(
    diamond_path,
    makedb_args,
    error_on_status = FALSE
  )

  if (makedb_result$status != 0) {
    cli::cli_abort("DIAMOND makedb failed: {makedb_result$stderr}")
  }

  # Run blastp
  blastp_args <- c(
    "blastp",
    "-d", db_path,
    "-q", fasta_path,
    "-o", hits_path,
    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore",
    "-p", as.character(threads),
    "-e", as.character(evalue),
    "--max-target-seqs", if (mode == "fast") "50" else "100",
    "--quiet"
  )

  # Add sensitivity for thorough mode
  if (mode == "thorough") {
    blastp_args <- c(blastp_args, "--sensitive")
  }

  blastp_result <- processx::run(
    diamond_path,
    blastp_args,
    error_on_status = FALSE
  )

  if (blastp_result$status != 0) {
    cli::cli_abort("DIAMOND blastp failed: {blastp_result$stderr}")
  }

  # Parse hits with filtering
  hits <- parse_diamond_hits(hits_path, min_identity, min_coverage, evalue)

  # Compute RBH
  rbh_pairs <- compute_rbh(hits)

  # Get all protein IDs
  all_protein_ids <- unlist(lapply(proteins$assemblies, function(ps) {
    ps$proteins$protein_id
  }))

  # Build orthogroups
  orthogroups_raw <- build_orthogroups_from_rbh(rbh_pairs, all_protein_ids)

  # Add assembly information to orthogroups
  # Create protein_id -> assembly lookup
  protein_assembly <- lapply(proteins$assemblies, function(ps) {
    tibble::tibble(
      protein_id = ps$proteins$protein_id,
      assembly = ps$assembly_name
    )
  }) |> dplyr::bind_rows()

  # Join assembly info
  orthogroups <- if (nrow(orthogroups_raw) > 0) {
    orthogroups_raw |>
      dplyr::left_join(protein_assembly, by = "protein_id") |>
      dplyr::select("orthogroup_id", "assembly", "protein_id")
  } else {
    tibble::tibble(
      orthogroup_id = character(),
      assembly = character(),
      protein_id = character()
    )
  }

  # Identify singletons
  clustered_proteins <- orthogroups$protein_id
  singleton_ids <- setdiff(all_protein_ids, clustered_proteins)

  singletons <- if (length(singleton_ids) > 0) {
    protein_assembly |>
      dplyr::filter(.data$protein_id %in% singleton_ids)
  } else {
    tibble::tibble(
      protein_id = character(),
      assembly = character()
    )
  }

  # Build result
  new_orthogroup_result(
    orthogroups = orthogroups,
    method = "diamond_rbh",
    parameters = list(
      min_identity = min_identity,
      min_coverage = min_coverage,
      evalue = evalue,
      mode = mode,
      threads = threads
    ),
    singletons = singletons
  )
}

# Main Dispatcher Function --------------------------------------------------

#' Cluster proteins across assemblies
#'
#' Groups proteins from a protein_collection into orthogroups using
#' sequence similarity. Multiple clustering methods are available.
#'
#' @param proteins A protein_collection object
#' @param method Clustering method: "diamond_rbh" (default), "orthofinder", or "mmseqs2"
#' @param mode Speed/sensitivity trade-off: "fast" (default) or "thorough"
#' @param min_identity Minimum percent identity threshold (default 30)
#' @param min_coverage Minimum query coverage threshold (default 50)
#' @param evalue E-value threshold (default 1e-5)
#' @param threads Number of CPU threads (default: auto-detect)
#' @param tool_path Optional explicit path to the clustering tool binary
#' @param conda_env Optional path to conda/mamba environment containing the tool
#' @param keep_temp Keep temporary files for debugging (default FALSE)
#'
#' @return An orthogroup_result object containing:
#'   - orthogroups: tibble with orthogroup_id, assembly, protein_id
#'   - method: the clustering method used
#'   - parameters: list of parameters used
#'   - singletons: proteins not assigned to any orthogroup
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load proteins
#' proteins <- load_proteins(fasta_dir = "assemblies/")
#'
#' # Cluster with default settings (DIAMOND RBH)
#' result <- cluster_proteins(proteins)
#'
#' # Use thorough mode with stricter thresholds
#' result <- cluster_proteins(
#'   proteins,
#'   mode = "thorough",
#'   min_identity = 70,
#'   min_coverage = 80
#' )
#'
#' # Use a conda environment
#' result <- cluster_proteins(
#'   proteins,
#'   conda_env = "./this_project_env"
#' )
#' }
cluster_proteins <- function(proteins,
                             method = "diamond_rbh",
                             mode = "fast",
                             min_identity = 30,
                             min_coverage = 50,
                             evalue = 1e-5,
                             threads = NULL,
                             tool_path = NULL,
                             conda_env = NULL,
                             keep_temp = FALSE) {
  # Validate proteins argument
 if (!inherits(proteins, "protein_collection")) {
    cli::cli_abort("{.arg proteins} must be a {.cls protein_collection} object")
  }

  # Validate method
  valid_methods <- c("diamond_rbh", "orthofinder", "mmseqs2")
  if (!method %in% valid_methods) {
    cli::cli_abort(
      "{.arg method} must be one of: {.val {valid_methods}}"
    )
  }

  # Validate mode
  valid_modes <- c("fast", "thorough")
  if (!mode %in% valid_modes) {
    cli::cli_abort(
      "{.arg mode} must be one of: {.val {valid_modes}}"
    )
  }

  # Build tool_path from conda_env if provided
  if (!is.null(conda_env) && is.null(tool_path)) {
    tool_path <- file.path(conda_env, "bin", "diamond")
  }

  # Dispatch to appropriate method
  if (method == "diamond_rbh") {
    run_diamond_rbh(
      proteins = proteins,
      min_identity = min_identity,
      min_coverage = min_coverage,
      evalue = evalue,
      threads = threads,
      mode = mode,
      tool_path = tool_path,
      keep_temp = keep_temp
    )
  } else if (method == "orthofinder") {
    cli::cli_abort("Method {.val orthofinder} is not yet implemented")
  } else if (method == "mmseqs2") {
    cli::cli_abort("Method {.val mmseqs2} is not yet implemented")
  }
}
