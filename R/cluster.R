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
