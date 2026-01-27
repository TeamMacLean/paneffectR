# Tests for S3 class constructors

# protein_set tests --------------------------------------------------------

test_that("new_protein_set creates valid object", {
  proteins <- tibble::tibble(
    protein_id = c("prot1", "prot2"),
    sequence = c("MAAAA", "MVVVV")
  )
  ps <- new_protein_set("assembly1", proteins)
  expect_s3_class(ps, "protein_set")
})

test_that("new_protein_set stores all components", {
  proteins <- tibble::tibble(
    protein_id = c("prot1", "prot2"),
    sequence = c("MAAAA", "MVVVV")
  )
  ps <- new_protein_set(
    "assembly1",
    proteins,
    metadata = list(species = "test"),
    source_files = list(fasta = "/path/to/file.faa")
  )
  expect_equal(ps$assembly_name, "assembly1")
  expect_equal(ps$proteins, proteins)
  expect_equal(ps$metadata$species, "test")
  expect_equal(ps$source_files$fasta, "/path/to/file.faa")
})

test_that("new_protein_set validates assembly_name", {
  proteins <- tibble::tibble(protein_id = "p1", sequence = "M")
  expect_error(new_protein_set(NULL, proteins))
  expect_error(new_protein_set("", proteins))
  expect_error(new_protein_set(123, proteins))
  expect_error(new_protein_set(c("a", "b"), proteins))
})

test_that("new_protein_set validates proteins tibble", {
  expect_error(new_protein_set("a", NULL))
  expect_error(new_protein_set("a", data.frame(x = 1)))
  expect_error(new_protein_set("a", tibble::tibble(x = 1)))
})

test_that("new_protein_set rejects duplicate protein_ids", {
  proteins <- tibble::tibble(
    protein_id = c("prot1", "prot1"),
    sequence = c("MAAAA", "MVVVV")
  )
  expect_error(new_protein_set("assembly1", proteins))
})

test_that("print.protein_set produces output", {
  proteins <- tibble::tibble(protein_id = "p1", sequence = "M")
  ps <- new_protein_set("assembly1", proteins)
  expect_output(print(ps), "assembly1")
  expect_output(print(ps), "1 protein")
})

# protein_collection tests -------------------------------------------------

test_that("new_protein_collection creates valid object", {
  ps1 <- new_protein_set("a1", tibble::tibble(protein_id = "p1", sequence = "M"))
  ps2 <- new_protein_set("a2", tibble::tibble(protein_id = "p2", sequence = "M"))
  pc <- new_protein_collection(list(ps1, ps2))
  expect_s3_class(pc, "protein_collection")
})

test_that("new_protein_collection computes n_assemblies", {
  ps1 <- new_protein_set("a1", tibble::tibble(protein_id = "p1", sequence = "M"))
  ps2 <- new_protein_set("a2", tibble::tibble(protein_id = "p2", sequence = "M"))
  pc <- new_protein_collection(list(ps1, ps2))
  expect_equal(pc$n_assemblies, 2)
})

test_that("new_protein_collection computes summary", {
  ps1 <- new_protein_set("a1", tibble::tibble(protein_id = "p1", sequence = "M"))
  ps2 <- new_protein_set("a2", tibble::tibble(protein_id = c("p2", "p3"), sequence = c("M", "M")))
  pc <- new_protein_collection(list(ps1, ps2))
  expect_s3_class(pc$summary, "tbl_df")
  expect_equal(nrow(pc$summary), 2)
  expect_true("assembly_name" %in% names(pc$summary))
  expect_true("n_proteins" %in% names(pc$summary))
})

test_that("new_protein_collection rejects duplicate assembly names", {
  ps1 <- new_protein_set("a1", tibble::tibble(protein_id = "p1", sequence = "M"))
  ps2 <- new_protein_set("a1", tibble::tibble(protein_id = "p2", sequence = "M"))
  expect_error(new_protein_collection(list(ps1, ps2)))
})

test_that("new_protein_collection rejects duplicate protein_ids across assemblies", {
  ps1 <- new_protein_set("a1", tibble::tibble(protein_id = "p1", sequence = "M"))
  ps2 <- new_protein_set("a2", tibble::tibble(protein_id = "p1", sequence = "M"))
  expect_error(new_protein_collection(list(ps1, ps2)))
})

test_that("new_protein_collection validates input is list of protein_sets", {
  expect_error(new_protein_collection(NULL))
  expect_error(new_protein_collection(list()))
  expect_error(new_protein_collection(list("not a protein_set")))
})

test_that("print.protein_collection produces output", {
  ps1 <- new_protein_set("a1", tibble::tibble(protein_id = "p1", sequence = "M"))
  pc <- new_protein_collection(list(ps1))
  expect_output(print(pc), "1 assembl")
})

# orthogroup_result and pa_matrix stubs ------------------------------------

test_that("new_orthogroup_result is not yet implemented", {
  expect_error(new_orthogroup_result(NULL, "diamond_rbh"), "not yet implemented")
})

test_that("new_pa_matrix is not yet implemented", {
  expect_error(new_pa_matrix(NULL, NULL, NULL), "not yet implemented")
})
