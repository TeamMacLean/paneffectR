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

# orthogroup_result tests ---------------------------------------------------

test_that("new_orthogroup_result creates valid object", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002", "OG0002"),
    assembly = c("a1", "a2", "a1", "a2"),
    protein_id = c("a1_p1", "a2_p1", "a1_p2", "a2_p2")
  )

  result <- new_orthogroup_result(orthogroups, method = "diamond_rbh")
  expect_s3_class(result, "orthogroup_result")
})

test_that("new_orthogroup_result stores all components", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001"),
    assembly = c("a1", "a2"),
    protein_id = c("a1_p1", "a2_p1")
  )
  singletons <- tibble::tibble(
    assembly = "a1",
    protein_id = "a1_p99"
  )
  params <- list(min_identity = 70, min_coverage = 50)

  result <- new_orthogroup_result(
    orthogroups,
    method = "diamond_rbh",
    parameters = params,
    singletons = singletons
  )

  expect_equal(result$orthogroups, orthogroups)
  expect_equal(result$method, "diamond_rbh")
  expect_equal(result$parameters, params)
  expect_equal(result$singletons, singletons)
})

test_that("new_orthogroup_result validates orthogroups is tibble", {
  expect_error(new_orthogroup_result(NULL, "test"))
  expect_error(new_orthogroup_result(data.frame(x = 1), "test"))
  expect_error(new_orthogroup_result("not a tibble", "test"))
})
test_that("new_orthogroup_result validates required columns", {
  # Missing orthogroup_id
  bad1 <- tibble::tibble(assembly = "a1", protein_id = "p1")
  expect_error(new_orthogroup_result(bad1, "test"), "orthogroup_id")

  # Missing assembly
  bad2 <- tibble::tibble(orthogroup_id = "OG1", protein_id = "p1")
  expect_error(new_orthogroup_result(bad2, "test"), "assembly")

  # Missing protein_id
  bad3 <- tibble::tibble(orthogroup_id = "OG1", assembly = "a1")
  expect_error(new_orthogroup_result(bad3, "test"), "protein_id")
})

test_that("new_orthogroup_result validates method is character", {
  orthogroups <- tibble::tibble(
    orthogroup_id = "OG1",
    assembly = "a1",
    protein_id = "p1"
  )
  expect_error(new_orthogroup_result(orthogroups, method = NULL))
  expect_error(new_orthogroup_result(orthogroups, method = 123))
  expect_error(new_orthogroup_result(orthogroups, method = c("a", "b")))
})

test_that("new_orthogroup_result validates parameters is list", {
  orthogroups <- tibble::tibble(
    orthogroup_id = "OG1",
    assembly = "a1",
    protein_id = "p1"
  )
  expect_error(new_orthogroup_result(orthogroups, "test", parameters = "not a list"))
})

test_that("new_orthogroup_result validates singletons if provided", {
  orthogroups <- tibble::tibble(
    orthogroup_id = "OG1",
    assembly = "a1",
    protein_id = "p1"
  )
  # Not a tibble
  expect_error(new_orthogroup_result(orthogroups, "test", singletons = "bad"))

  # Missing required columns
  bad_singletons <- tibble::tibble(x = 1)
  expect_error(new_orthogroup_result(orthogroups, "test", singletons = bad_singletons))
})

test_that("new_orthogroup_result computes stats automatically", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001", "OG0002", "OG0002"),
    assembly = c("a1", "a2", "a3", "a1", "a2"),
    protein_id = c("a1_p1", "a2_p1", "a3_p1", "a1_p2", "a2_p2")
  )
  singletons <- tibble::tibble(
    assembly = c("a3", "a3"),
    protein_id = c("a3_p2", "a3_p3")
  )

  result <- new_orthogroup_result(orthogroups, "test", singletons = singletons)

  expect_true("stats" %in% names(result))
  expect_s3_class(result$stats, "tbl_df")

  # Should have n_orthogroups
  expect_true("n_orthogroups" %in% names(result$stats))
  expect_equal(result$stats$n_orthogroups, 2)

  # Should have n_singletons
  expect_true("n_singletons" %in% names(result$stats))
  expect_equal(result$stats$n_singletons, 2)

  # Should have n_proteins_clustered
  expect_true("n_proteins_clustered" %in% names(result$stats))
  expect_equal(result$stats$n_proteins_clustered, 5)

  # Should have n_assemblies
  expect_true("n_assemblies" %in% names(result$stats))
  expect_equal(result$stats$n_assemblies, 3)
})

test_that("new_orthogroup_result stats handles no singletons", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG1", "OG1"),
    assembly = c("a1", "a2"),
    protein_id = c("p1", "p2")
  )
  result <- new_orthogroup_result(orthogroups, "test")
  expect_equal(result$stats$n_singletons, 0)
})

test_that("new_orthogroup_result allows user-provided stats", {
 orthogroups <- tibble::tibble(
    orthogroup_id = "OG1",
    assembly = "a1",
    protein_id = "p1"
  )
  custom_stats <- tibble::tibble(
    n_orthogroups = 1L,
    n_singletons = 0L,
    n_proteins_clustered = 1L,
    n_assemblies = 1L,
    custom_field = "extra"
  )
  result <- new_orthogroup_result(orthogroups, "test", stats = custom_stats)
  expect_equal(result$stats$custom_field, "extra")
})

test_that("print.orthogroup_result produces output", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002"),
    assembly = c("a1", "a2", "a1"),
    protein_id = c("p1", "p2", "p3")
  )
  result <- new_orthogroup_result(orthogroups, "diamond_rbh")

  expect_output(print(result), "orthogroup_result")
  expect_output(print(result), "diamond_rbh")
  expect_output(print(result), "2 orthogroup")
})

test_that("print.orthogroup_result shows singletons when present", {
  orthogroups <- tibble::tibble(
    orthogroup_id = "OG1",
    assembly = "a1",
    protein_id = "p1"
  )
  singletons <- tibble::tibble(
    assembly = c("a1", "a2"),
    protein_id = c("s1", "s2")
  )
  result <- new_orthogroup_result(orthogroups, "test", singletons = singletons)
  expect_output(print(result), "2 singleton")
})

# pa_matrix tests -----------------------------------------------------------

test_that("new_pa_matrix creates valid object", {
  mat <- matrix(
    c(1, 0, 1, 1, 1, 0),
    nrow = 3, ncol = 2,
    dimnames = list(c("OG0001", "OG0002", "OG0003"), c("asm1", "asm2"))
  )
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002", "OG0003"),
    size = c(2L, 2L, 1L)
  )
  assemblies <- tibble::tibble(
    assembly_name = c("asm1", "asm2"),
    n_orthogroups = c(3L, 2L)
  )

  pa <- new_pa_matrix(mat, orthogroups, assemblies, type = "binary")
  expect_s3_class(pa, "pa_matrix")
})

test_that("new_pa_matrix stores all components", {
  mat <- matrix(
    c(1, 0, 1, 1),
    nrow = 2, ncol = 2,
    dimnames = list(c("OG0001", "OG0002"), c("asm1", "asm2"))
  )
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    size = c(2L, 2L)
  )
  assemblies <- tibble::tibble(
    assembly_name = c("asm1", "asm2"),
    n_orthogroups = c(2L, 2L)
  )

  pa <- new_pa_matrix(mat, orthogroups, assemblies, type = "binary", threshold = 0.5)

  expect_equal(pa$matrix, mat)
  expect_equal(pa$orthogroups, orthogroups)
  expect_equal(pa$assemblies, assemblies)
  expect_equal(pa$type, "binary")
  expect_equal(pa$threshold, 0.5)
})

test_that("new_pa_matrix validates matrix is a matrix", {
  orthogroups <- tibble::tibble(orthogroup_id = "OG1", size = 1L)
  assemblies <- tibble::tibble(assembly_name = "asm1", n_orthogroups = 1L)

  expect_error(new_pa_matrix(NULL, orthogroups, assemblies, "binary"))
  expect_error(new_pa_matrix("not a matrix", orthogroups, assemblies, "binary"))
  expect_error(new_pa_matrix(data.frame(x = 1), orthogroups, assemblies, "binary"))
})

test_that("new_pa_matrix validates matrix row count matches orthogroups", {
  # Matrix has 3 rows, orthogroups has 2
  mat <- matrix(
    c(1, 0, 1, 1, 1, 0),
    nrow = 3, ncol = 2,
    dimnames = list(c("OG0001", "OG0002", "OG0003"), c("asm1", "asm2"))
  )
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    size = c(2L, 2L)
  )
  assemblies <- tibble::tibble(
    assembly_name = c("asm1", "asm2"),
    n_orthogroups = c(2L, 2L)
  )

  expect_error(new_pa_matrix(mat, orthogroups, assemblies, "binary"), "row")
})

test_that("new_pa_matrix validates matrix column count matches assemblies", {
  # Matrix has 2 cols, assemblies has 3
  mat <- matrix(
    c(1, 0, 1, 1),
    nrow = 2, ncol = 2,
    dimnames = list(c("OG0001", "OG0002"), c("asm1", "asm2"))
  )
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    size = c(2L, 2L)
  )
  assemblies <- tibble::tibble(
    assembly_name = c("asm1", "asm2", "asm3"),
    n_orthogroups = c(2L, 2L, 1L)
  )

  expect_error(new_pa_matrix(mat, orthogroups, assemblies, "binary"), "column")
})

test_that("new_pa_matrix validates orthogroups is tibble with required columns", {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("OG1", "asm1"))
  assemblies <- tibble::tibble(assembly_name = "asm1", n_orthogroups = 1L)

  # Not a tibble
  expect_error(new_pa_matrix(mat, NULL, assemblies, "binary"))
  expect_error(new_pa_matrix(mat, "not tibble", assemblies, "binary"))

  # Missing orthogroup_id column
  bad_og <- tibble::tibble(size = 1L)
  expect_error(new_pa_matrix(mat, bad_og, assemblies, "binary"), "orthogroup_id")
})

test_that("new_pa_matrix validates assemblies is tibble with required columns", {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("OG1", "asm1"))
  orthogroups <- tibble::tibble(orthogroup_id = "OG1", size = 1L)

  # Not a tibble
  expect_error(new_pa_matrix(mat, orthogroups, NULL, "binary"))
  expect_error(new_pa_matrix(mat, orthogroups, "not tibble", "binary"))

  # Missing assembly_name column
  bad_asm <- tibble::tibble(n_orthogroups = 1L)
  expect_error(new_pa_matrix(mat, orthogroups, bad_asm, "binary"), "assembly_name")
})

test_that("new_pa_matrix validates type is one of allowed values", {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("OG1", "asm1"))
  orthogroups <- tibble::tibble(orthogroup_id = "OG1", size = 1L)
  assemblies <- tibble::tibble(assembly_name = "asm1", n_orthogroups = 1L)

  expect_error(new_pa_matrix(mat, orthogroups, assemblies, "invalid"))
  expect_error(new_pa_matrix(mat, orthogroups, assemblies, NULL))
  expect_error(new_pa_matrix(mat, orthogroups, assemblies, 123))

  # Valid types should work (once implemented)
  # expect_no_error(new_pa_matrix(mat, orthogroups, assemblies, "binary"))
  # expect_no_error(new_pa_matrix(mat, orthogroups, assemblies, "count"))
  # expect_no_error(new_pa_matrix(mat, orthogroups, assemblies, "score"))
})

test_that("new_pa_matrix accepts NULL threshold for binary and count types", {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("OG1", "asm1"))
  orthogroups <- tibble::tibble(orthogroup_id = "OG1", size = 1L)
  assemblies <- tibble::tibble(assembly_name = "asm1", n_orthogroups = 1L)

  pa_binary <- new_pa_matrix(mat, orthogroups, assemblies, "binary", threshold = NULL)
  expect_null(pa_binary$threshold)

  pa_count <- new_pa_matrix(mat, orthogroups, assemblies, "count", threshold = NULL)
  expect_null(pa_count$threshold)
})

test_that("new_pa_matrix accepts numeric threshold", {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("OG1", "asm1"))
  orthogroups <- tibble::tibble(orthogroup_id = "OG1", size = 1L)
  assemblies <- tibble::tibble(assembly_name = "asm1", n_orthogroups = 1L)

  pa <- new_pa_matrix(mat, orthogroups, assemblies, "score", threshold = 5.0)
  expect_equal(pa$threshold, 5.0)
})

test_that("print.pa_matrix shows dimensions", {
  mat <- matrix(
    c(1, 0, 1, 1, 1, 0),
    nrow = 3, ncol = 2,
    dimnames = list(c("OG0001", "OG0002", "OG0003"), c("asm1", "asm2"))
  )
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002", "OG0003"),
    size = c(2L, 2L, 1L)
  )
  assemblies <- tibble::tibble(
    assembly_name = c("asm1", "asm2"),
    n_orthogroups = c(3L, 2L)
  )
  pa <- new_pa_matrix(mat, orthogroups, assemblies, "binary")

  expect_output(print(pa), "3")  # 3 orthogroups
  expect_output(print(pa), "2")  # 2 assemblies
})

test_that("print.pa_matrix shows type", {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("OG1", "asm1"))
  orthogroups <- tibble::tibble(orthogroup_id = "OG1", size = 1L)
  assemblies <- tibble::tibble(assembly_name = "asm1", n_orthogroups = 1L)
  pa <- new_pa_matrix(mat, orthogroups, assemblies, "binary")

  expect_output(print(pa), "binary")
})

test_that("print.pa_matrix shows sparsity percentage", {
  # Matrix with 2 zeros out of 6 cells = 33.3% sparse
  mat <- matrix(
    c(1, 0, 1, 1, 1, 0),
    nrow = 3, ncol = 2,
    dimnames = list(c("OG0001", "OG0002", "OG0003"), c("asm1", "asm2"))
  )
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002", "OG0003"),
    size = c(2L, 2L, 1L)
  )
  assemblies <- tibble::tibble(
    assembly_name = c("asm1", "asm2"),
    n_orthogroups = c(3L, 2L)
  )
  pa <- new_pa_matrix(mat, orthogroups, assemblies, "binary")

  # Should show some percentage indicating sparsity
  expect_output(print(pa), "%")
})
