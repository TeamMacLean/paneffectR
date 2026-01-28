# Tests for clustering functions (Phase 2)

# compute_rbh tests ---------------------------------------------------------

test_that("compute_rbh finds reciprocal best hits", {
 # A_001's best hit is B_001 (score 100)
  # B_001's best hit is A_001 (score 100)
  # This is a reciprocal best hit pair
  hits <- tibble::tibble(
    qseqid = c("A_001", "A_001", "B_001", "B_002"),
    sseqid = c("B_001", "B_002", "A_001", "A_001"),
    bitscore = c(100, 50, 100, 80)
  )
  rbh <- compute_rbh(hits)

  expect_s3_class(rbh, "tbl_df")
  expect_equal(nrow(rbh), 1)
  expect_true(all(c("protein_a", "protein_b") %in% names(rbh)))
  # The pair should be A_001 and B_001
 expect_true(
    (rbh$protein_a[1] == "A_001" && rbh$protein_b[1] == "B_001") ||
    (rbh$protein_a[1] == "B_001" && rbh$protein_b[1] == "A_001")
  )
})

test_that("compute_rbh handles multiple RBH pairs", {
  # Two reciprocal pairs: A_001<->B_001 and A_002<->B_002
  hits <- tibble::tibble(
    qseqid = c("A_001", "A_002", "B_001", "B_002"),
    sseqid = c("B_001", "B_002", "A_001", "A_002"),
    bitscore = c(100, 90, 100, 90)
  )
  rbh <- compute_rbh(hits)
  expect_equal(nrow(rbh), 2)
})

test_that("compute_rbh excludes non-reciprocal hits", {
  # A_001's best is B_001, but B_001's best is A_002 (not reciprocal)
  hits <- tibble::tibble(
    qseqid = c("A_001", "A_002", "B_001"),
    sseqid = c("B_001", "B_001", "A_002"),
    bitscore = c(80, 100, 100)
  )
  rbh <- compute_rbh(hits)
  # Only A_002<->B_001 is reciprocal
  expect_equal(nrow(rbh), 1)
  expect_true("A_002" %in% c(rbh$protein_a[1], rbh$protein_b[1]))
})

test_that("compute_rbh excludes self-hits", {
  hits <- tibble::tibble(
    qseqid = c("A_001", "A_001", "B_001"),
    sseqid = c("A_001", "B_001", "A_001"),
    bitscore = c(200, 100, 100)
  )
  rbh <- compute_rbh(hits)
  # Self-hit should be excluded, only A_001<->B_001 remains
  expect_equal(nrow(rbh), 1)
})

test_that("compute_rbh returns empty tibble when no RBH found", {
  # No reciprocal relationships
  hits <- tibble::tibble(
    qseqid = c("A_001", "B_001"),
    sseqid = c("B_001", "C_001"),
    bitscore = c(100, 100)
  )
  rbh <- compute_rbh(hits)
  expect_equal(nrow(rbh), 0)
  expect_true(all(c("protein_a", "protein_b") %in% names(rbh)))
})

test_that("compute_rbh handles empty input", {
  hits <- tibble::tibble(
    qseqid = character(),
    sseqid = character(),
    bitscore = numeric()
  )
  rbh <- compute_rbh(hits)
  expect_equal(nrow(rbh), 0)
})

# build_orthogroups_from_rbh tests ------------------------------------------

test_that("build_orthogroups_from_rbh creates connected components", {
  # A-B and B-C should form one group (A, B, C connected)
  rbh <- tibble::tibble(
    protein_a = c("A", "B"),
    protein_b = c("B", "C")
  )
  result <- build_orthogroups_from_rbh(rbh, c("A", "B", "C", "D"))

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("orthogroup_id", "protein_id") %in% names(result)))

  # A, B, C should be in the same orthogroup
  og_for_abc <- result$orthogroup_id[result$protein_id %in% c("A", "B", "C")]
  expect_equal(length(unique(og_for_abc)), 1)

  # D is a singleton, should not be in result
  expect_false("D" %in% result$protein_id)
})

test_that("build_orthogroups_from_rbh handles multiple components", {
  # Two separate components: A-B and C-D
  rbh <- tibble::tibble(
    protein_a = c("A", "C"),
    protein_b = c("B", "D")
  )
  result <- build_orthogroups_from_rbh(rbh, c("A", "B", "C", "D", "E"))

  # Should have 2 orthogroups
  expect_equal(length(unique(result$orthogroup_id)), 2)

  # A and B in same group
  expect_equal(
    result$orthogroup_id[result$protein_id == "A"],
    result$orthogroup_id[result$protein_id == "B"]
  )

  # C and D in same group
  expect_equal(
    result$orthogroup_id[result$protein_id == "C"],
    result$orthogroup_id[result$protein_id == "D"]
  )

  # E is singleton, not in result
  expect_false("E" %in% result$protein_id)
})

test_that("build_orthogroups_from_rbh returns singletons separately",
{
  rbh <- tibble::tibble(
    protein_a = "A",
    protein_b = "B"
  )
  all_proteins <- c("A", "B", "C", "D")

  result <- build_orthogroups_from_rbh(rbh, all_proteins)

  # Result should only contain clustered proteins
  expect_equal(sort(result$protein_id), c("A", "B"))
})

test_that("build_orthogroups_from_rbh handles empty RBH", {
  rbh <- tibble::tibble(
    protein_a = character(),
    protein_b = character()
  )
  result <- build_orthogroups_from_rbh(rbh, c("A", "B", "C"))

  # No orthogroups when no RBH pairs
  expect_equal(nrow(result), 0)
})

test_that("build_orthogroups_from_rbh generates sequential OG IDs", {
  rbh <- tibble::tibble(
    protein_a = c("A", "C"),
    protein_b = c("B", "D")
  )
  result <- build_orthogroups_from_rbh(rbh, c("A", "B", "C", "D"))

  # IDs should be like OG0001, OG0002, etc.
  expect_true(all(grepl("^OG\\d+$", result$orthogroup_id)))
})
