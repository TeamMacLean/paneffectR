# Tests for data loading functions

# Helper to get test data path
get_testdata_path <- function(file = NULL) {
  path <- system.file("testdata", package = "paneffectR")
  if (path == "") {
    # Fallback for development when package isn't installed
    path <- test_path("../../inst/testdata")
  }
  if (!is.null(file)) {
    path <- file.path(path, file)
  }
  path
}

# load_fasta() tests -------------------------------------------------------

test_that("load_fasta returns tibble with protein_id and sequence columns", {
  path <- get_testdata_path("assembly1.faa")
  result <- load_fasta(path)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("protein_id", "sequence"))
})

test_that("load_fasta parses correct number of proteins", {
  path <- get_testdata_path("assembly1.faa")
  result <- load_fasta(path)
  expect_equal(nrow(result), 100)
})

test_that("load_fasta extracts protein IDs correctly", {
  path <- get_testdata_path("assembly1.faa")
  result <- load_fasta(path)
  # Check some known protein IDs from the test data
  expect_true(all(grepl("^assembly1_", result$protein_id)))
})

test_that("load_fasta handles multi-line sequences", {
  path <- get_testdata_path("assembly1.faa")
  result <- load_fasta(path)
  # Sequences should be concatenated, not split across lines
  # All sequences should be non-empty strings without newlines
  expect_true(all(nchar(result$sequence) > 0))
  expect_false(any(grepl("\n", result$sequence)))
})

test_that("load_fasta errors on non-existent file", {
  expect_error(load_fasta("/nonexistent/path.faa"))
})

test_that("load_fasta errors on invalid path type", {
  expect_error(load_fasta(123))
  expect_error(load_fasta(NULL))
})

# load_scores() tests ------------------------------------------------------

test_that("load_scores returns tibble with expected columns", {
  path <- get_testdata_path("assembly1_scored.csv")
  result <- load_scores(path)
  expect_s3_class(result, "tbl_df")
  expect_true("assembly" %in% names(result))
  expect_true("protein_id" %in% names(result))
  expect_true("custom_score" %in% names(result))
  expect_true("score_rank" %in% names(result))
})

test_that("load_scores parses correct number of rows", {
  path <- get_testdata_path("assembly1_scored.csv")
  result <- load_scores(path)
  expect_equal(nrow(result), 100)
})

test_that("load_scores parses numeric columns correctly", {
  path <- get_testdata_path("assembly1_scored.csv")
  result <- load_scores(path)
  expect_type(result$custom_score, "double")
  expect_type(result$score_rank, "integer")
})

test_that("load_scores errors on non-existent file", {
  expect_error(load_scores("/nonexistent/path.csv"))
})

# load_name_mapping() tests ------------------------------------------------

test_that("load_name_mapping returns tibble with required columns", {
  path <- get_testdata_path("assembly1_name_mapping.csv")
  result <- load_name_mapping(path)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("assembly", "old_name", "new_name", "protein_number") %in% names(result)))
})

test_that("load_name_mapping parses correct number of rows", {
  path <- get_testdata_path("assembly1_name_mapping.csv")
  result <- load_name_mapping(path)
  expect_equal(nrow(result), 100)
})

test_that("load_name_mapping errors on non-existent file", {
  expect_error(load_name_mapping("/nonexistent/path.csv"))
})

# load_proteins() tests ----------------------------------------------------

test_that("load_proteins loads single assembly from directory", {
  fasta_dir <- get_testdata_path()
  pc <- load_proteins(fasta_dir, pattern = "assembly1.faa")
  expect_s3_class(pc, "protein_collection")
  expect_equal(pc$n_assemblies, 1)
})

test_that("load_proteins loads multiple assemblies", {
  fasta_dir <- get_testdata_path()
  pc <- load_proteins(fasta_dir, pattern = "*.faa")
  expect_s3_class(pc, "protein_collection")
  expect_equal(pc$n_assemblies, 3)
})

test_that("load_proteins merges scores when score_dir provided", {
  fasta_dir <- get_testdata_path()
  score_dir <- get_testdata_path()
  pc <- load_proteins(fasta_dir, score_dir = score_dir, pattern = "assembly1.faa")
  # Check that score columns are present in proteins
  expect_true("custom_score" %in% names(pc$assemblies[[1]]$proteins))
})

test_that("load_proteins extracts assembly name from filename", {
  fasta_dir <- get_testdata_path()
  pc <- load_proteins(fasta_dir, pattern = "assembly1.faa")
  expect_equal(pc$assemblies[[1]]$assembly_name, "assembly1")
})

test_that("load_proteins errors on non-existent directory", {
  expect_error(load_proteins("/nonexistent/directory"))
})

test_that("load_proteins errors when no FASTA files match pattern", {
  fasta_dir <- get_testdata_path()
  expect_error(load_proteins(fasta_dir, pattern = "*.xyz"))
})
