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

# write_combined_fasta tests ------------------------------------------------
# TODO: Currently uses protein_id as header. May need assembly prefix
# if users have duplicate protein_ids across assemblies.

test_that("write_combined_fasta writes all proteins", {
  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = c("asm1_001", "asm1_002"),
    sequence = c("MAAAA", "MVVVV")
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = c("asm2_001"),
    sequence = c("MGGGG")
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  tmp <- tempfile(fileext = ".faa")
  on.exit(unlink(tmp))

  write_combined_fasta(pc, tmp)

  expect_true(file.exists(tmp))
  lines <- readLines(tmp)

  # Should have 3 proteins (6 lines: 3 headers + 3 sequences)
  headers <- grep("^>", lines, value = TRUE)
  expect_equal(length(headers), 3)

  # All protein IDs should be present in headers
  expect_true(any(grepl("asm1_001", headers)))
  expect_true(any(grepl("asm1_002", headers)))
  expect_true(any(grepl("asm2_001", headers)))
})

test_that("write_combined_fasta writes correct sequences", {
  ps <- new_protein_set("asm1", tibble::tibble(
    protein_id = c("p1", "p2"),
    sequence = c("MTEST", "MFOO")
  ))
  pc <- new_protein_collection(list(ps))

  tmp <- tempfile(fileext = ".faa")
  on.exit(unlink(tmp))

  write_combined_fasta(pc, tmp)

  lines <- readLines(tmp)
  expect_true("MTEST" %in% lines)
  expect_true("MFOO" %in% lines)
})

test_that("write_combined_fasta handles single assembly", {
  ps <- new_protein_set("asm1", tibble::tibble(
    protein_id = "p1",
    sequence = "MTEST"
  ))
  pc <- new_protein_collection(list(ps))

  tmp <- tempfile(fileext = ".faa")
  on.exit(unlink(tmp))

  write_combined_fasta(pc, tmp)

  lines <- readLines(tmp)
  expect_equal(length(grep("^>", lines)), 1)
  expect_true("MTEST" %in% lines)
})

# parse_diamond_hits tests --------------------------------------------------
# DIAMOND output format: qseqid sseqid pident length qlen evalue bitscore
# Coverage is computed as (length / qlen) * 100

test_that("parse_diamond_hits reads tabular format", {
  tmp <- tempfile()
  on.exit(unlink(tmp))

  # Format: qseqid sseqid pident length qlen evalue bitscore
  writeLines(c(
    "A_001\tB_001\t95.5\t90\t100\t1e-50\t200",
    "A_001\tB_002\t80.0\t85\t100\t1e-30\t150",
    "B_001\tA_001\t95.5\t90\t100\t1e-50\t200"
  ), tmp)

  hits <- parse_diamond_hits(tmp, min_identity = 70, min_coverage = 50, evalue = 1e-5)

  expect_s3_class(hits, "tbl_df")
  expect_true(all(c("qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore", "qcov") %in% names(hits)))
  expect_equal(nrow(hits), 3)
})

test_that("parse_diamond_hits computes coverage correctly", {
  tmp <- tempfile()
  on.exit(unlink(tmp))

  # length=80, qlen=100 -> coverage = 80%
  writeLines("A_001\tB_001\t95.0\t80\t100\t1e-50\t200", tmp)

  hits <- parse_diamond_hits(tmp, min_identity = 0, min_coverage = 0, evalue = 1)

  expect_equal(hits$qcov[1], 80)
})

test_that("parse_diamond_hits filters by identity", {
  tmp <- tempfile()
  on.exit(unlink(tmp))

  writeLines(c(
    "A_001\tB_001\t95.0\t90\t100\t1e-50\t200",
    "A_002\tB_002\t60.0\t90\t100\t1e-50\t150"
  ), tmp)

  hits <- parse_diamond_hits(tmp, min_identity = 70, min_coverage = 50, evalue = 1e-5)

  # Only the 95% identity hit should pass
  expect_equal(nrow(hits), 1)
  expect_equal(hits$qseqid[1], "A_001")
})

test_that("parse_diamond_hits filters by coverage", {
  tmp <- tempfile()
  on.exit(unlink(tmp))

  # First: length=90, qlen=100 -> 90% coverage (passes)
  # Second: length=30, qlen=100 -> 30% coverage (fails)
  writeLines(c(
    "A_001\tB_001\t95.0\t90\t100\t1e-50\t200",
    "A_002\tB_002\t95.0\t30\t100\t1e-50\t150"
  ), tmp)

  hits <- parse_diamond_hits(tmp, min_identity = 70, min_coverage = 50, evalue = 1e-5)

  # Only the 90% coverage hit should pass
  expect_equal(nrow(hits), 1)
  expect_equal(hits$qseqid[1], "A_001")
})

test_that("parse_diamond_hits filters by evalue", {
  tmp <- tempfile()
  on.exit(unlink(tmp))

  writeLines(c(
    "A_001\tB_001\t95.0\t90\t100\t1e-50\t200",
    "A_002\tB_002\t95.0\t90\t100\t0.1\t50"
  ), tmp)

  hits <- parse_diamond_hits(tmp, min_identity = 70, min_coverage = 50, evalue = 1e-5)

  # Only the 1e-50 evalue hit should pass
  expect_equal(nrow(hits), 1)
  expect_equal(hits$qseqid[1], "A_001")
})

test_that("parse_diamond_hits handles empty file", {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  file.create(tmp)

  hits <- parse_diamond_hits(tmp, min_identity = 70, min_coverage = 50, evalue = 1e-5)

  expect_s3_class(hits, "tbl_df")
  expect_equal(nrow(hits), 0)
  expect_true(all(c("qseqid", "sseqid", "pident", "qcov", "evalue", "bitscore") %in% names(hits)))
})

# run_diamond_rbh integration tests -----------------------------------------

# Helper to find diamond in project-local env (works from test directory)
get_diamond_path <- function() {
  # Try project-local path (from tests/testthat directory)
  local_path <- test_path("../../this_project_env/bin/diamond")
  if (file.exists(local_path)) return(local_path)

  # Try system PATH
  sys_path <- Sys.which("diamond")
  if (nzchar(sys_path)) return(sys_path)

  NULL
}

test_that("run_diamond_rbh returns orthogroup_result", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")

  # Create test protein_collection with some identical sequences
  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = c("asm1_001", "asm1_002"),
    sequence = c(
      "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY",
      "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD"
    )
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = c("asm2_001", "asm2_002"),
    sequence = c(
      # Same as asm1_001 - should cluster together
      "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY",
      "MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHG"
    )
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  result <- run_diamond_rbh(
    pc,
    min_identity = 90,
    min_coverage = 80,
    evalue = 1e-5,
    tool_path = diamond_path
  )

  expect_s3_class(result, "orthogroup_result")
  expect_equal(result$method, "diamond_rbh")
})

test_that("run_diamond_rbh clusters identical sequences", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")

  # Two assemblies with one identical protein each
  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = "asm1_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = "asm2_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  result <- run_diamond_rbh(
    pc,
    min_identity = 90,
    min_coverage = 80,
    evalue = 1e-5,
    tool_path = diamond_path
  )

  # Both proteins should be in the same orthogroup
  og_data <- result$orthogroups
  expect_equal(nrow(og_data), 2)
  expect_equal(length(unique(og_data$orthogroup_id)), 1)
})

test_that("run_diamond_rbh identifies singletons", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")

  # Two very different proteins - should not cluster
  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = "asm1_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = "asm2_001",
    sequence = "MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHG"
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  result <- run_diamond_rbh(
    pc,
    min_identity = 90,
    min_coverage = 80,
    evalue = 1e-5,
    tool_path = diamond_path
  )

  # Should have singletons (proteins not in any orthogroup)
  expect_true(!is.null(result$singletons))
  expect_equal(nrow(result$singletons), 2)
})

# cluster_proteins() dispatcher tests ----------------------------------------

test_that("cluster_proteins validates proteins argument", {
  # Not a protein_collection

  expect_error(
    cluster_proteins("not a collection"),
    "protein_collection"

  )

  expect_error(
    cluster_proteins(list(a = 1, b = 2)),
    "protein_collection"
  )

  expect_error(
    cluster_proteins(NULL),
    "protein_collection"
  )
})

test_that("cluster_proteins validates method argument", {
  ps <- new_protein_set("asm1", tibble::tibble(
    protein_id = "p1",
    sequence = "MTEST"
  ))
  pc <- new_protein_collection(list(ps))

  expect_error(
    cluster_proteins(pc, method = "invalid_method"),
    "method"
  )

  expect_error(
    cluster_proteins(pc, method = "blast"),
    "method"
  )
})
test_that("cluster_proteins validates mode argument", {
  ps <- new_protein_set("asm1", tibble::tibble(
    protein_id = "p1",
    sequence = "MTEST"
  ))
  pc <- new_protein_collection(list(ps))

  expect_error(
    cluster_proteins(pc, mode = "invalid_mode"),
    "mode"
  )

  expect_error(
    cluster_proteins(pc, mode = "ultra"),
    "mode"
  )
})

test_that("cluster_proteins errors for unimplemented methods", {
  ps <- new_protein_set("asm1", tibble::tibble(
    protein_id = "p1",
    sequence = "MTEST"
  ))
  pc <- new_protein_collection(list(ps))

  expect_error(
    cluster_proteins(pc, method = "orthofinder"),
    "not yet implemented"
  )

  expect_error(
    cluster_proteins(pc, method = "mmseqs2"),
    "not yet implemented"
  )
})

test_that("cluster_proteins dispatches to diamond_rbh", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")

  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = "asm1_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = "asm2_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  result <- cluster_proteins(
    pc,
    method = "diamond_rbh",
    min_identity = 90,
    min_coverage = 80,
    tool_path = diamond_path
  )

  expect_s3_class(result, "orthogroup_result")
  expect_equal(result$method, "diamond_rbh")
  # Both identical proteins should cluster
  expect_equal(nrow(result$orthogroups), 2)
  expect_equal(length(unique(result$orthogroups$orthogroup_id)), 1)
})

test_that("cluster_proteins uses conda_env to find tool", {
  # Skip if the project env doesn't exist
  env_path <- test_path("../../this_project_env")
  skip_if_not(dir.exists(env_path), "Project conda env not found")

  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = "asm1_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = "asm2_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  # Use conda_env instead of tool_path
  result <- cluster_proteins(
    pc,
    method = "diamond_rbh",
    min_identity = 90,
    min_coverage = 80,
    conda_env = env_path
  )

  expect_s3_class(result, "orthogroup_result")
  expect_equal(result$method, "diamond_rbh")
})

test_that("cluster_proteins passes mode parameter correctly", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")

  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = "asm1_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = "asm2_001",
    sequence = "MKTIIALSYIFCLVFADYKNTDNEANEPSDPYIKKATSALTKSYTNSDKGKLIKAATTTAKQSGKY"
  ))
  pc <- new_protein_collection(list(ps1, ps2))

  # Test thorough mode runs without error
  result <- cluster_proteins(
    pc,
    method = "diamond_rbh",
    mode = "thorough",
    min_identity = 90,
    min_coverage = 80,
    tool_path = diamond_path
  )

  expect_s3_class(result, "orthogroup_result")
  # Mode should be recorded in parameters
  expect_equal(result$parameters$mode, "thorough")
})
