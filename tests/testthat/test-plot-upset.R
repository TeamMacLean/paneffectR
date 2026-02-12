# Tests for plot_upset()
# Phase 4, Chunk 3: UpSet Plot
# Spec: specs/phase4-chunk3-upset-plot.md

# Helper to create test data
make_test_pa_matrix <- function() {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001",
                      "OG0002", "OG0002",
                      "OG0003"),
    assembly = c("asm1", "asm2", "asm3",
                 "asm1", "asm2",
                 "asm1"),
    protein_id = c("asm1_p1", "asm2_p1", "asm3_p1",
                   "asm1_p2", "asm2_p2",
                   "asm1_p3")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  build_pa_matrix(ort, exclude_singletons = TRUE)
}

# plot_upset() basic tests ----------------------------------------------------

test_that("plot_upset returns without error on valid pa_matrix", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  # Should not error
  expect_no_error(plot_upset(pa))
})

test_that("plot_upset validates input is pa_matrix", {
  skip_if_not_installed("UpSetR")

  expect_error(
    plot_upset("not a pa_matrix"),
    regexp = "pa_matrix"
  )
  expect_error(
    plot_upset(NULL),
    regexp = "pa_matrix"
  )
})

test_that("plot_upset works with binary matrix", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  expect_no_error(plot_upset(pa))
})

test_that("plot_upset works with count matrix", {
  skip_if_not_installed("UpSetR")

  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001"),
    assembly = c("asm1", "asm1", "asm2"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, type = "count", exclude_singletons = TRUE)

  # Should convert to binary internally
  expect_no_error(plot_upset(pa))
})

test_that("plot_upset works with score matrix", {
  skip_if_not_installed("UpSetR")

  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001"),
    assembly = c("asm1", "asm2"),
    protein_id = c("asm1_p1", "asm2_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  proteins_asm1 <- tibble::tibble(
    protein_id = "asm1_p1", sequence = "AAA", custom_score = 5
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = "asm2_p1", sequence = "BBB", custom_score = 7
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  proteins <- new_protein_collection(list(ps1, ps2))

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        exclude_singletons = TRUE)

  # Should convert to binary (NA -> 0, score -> 1)
  expect_no_error(plot_upset(pa))
})

# plot_upset() parameter tests ------------------------------------------------

test_that("plot_upset min_size parameter filters intersections", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  # Should not error with min_size
  expect_no_error(plot_upset(pa, min_size = 1))
  expect_no_error(plot_upset(pa, min_size = 2))
})

test_that("plot_upset max_sets parameter limits assemblies", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  # Should not error with max_sets
  expect_no_error(plot_upset(pa, max_sets = 2))
})

test_that("plot_upset order_by = 'freq' works", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  expect_no_error(plot_upset(pa, order_by = "freq"))
})

test_that("plot_upset order_by = 'degree' works", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  expect_no_error(plot_upset(pa, order_by = "degree"))
})

test_that("plot_upset errors on invalid order_by", {
  skip_if_not_installed("UpSetR")

  pa <- make_test_pa_matrix()

  expect_error(
    plot_upset(pa, order_by = "invalid"),
    regexp = "order_by"
  )
})

# plot_upset() edge cases -----------------------------------------------------

test_that("plot_upset handles single assembly with informative message", {
  skip_if_not_installed("UpSetR")

  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    assembly = c("asm1", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, exclude_singletons = TRUE)

  # Should error or warn - UpSet needs multiple sets
  expect_error(
    plot_upset(pa),
    regexp = "assembl|set"
  )
})

test_that("plot_upset errors on empty matrix", {
  skip_if_not_installed("UpSetR")

  # Create empty pa_matrix
  mat <- matrix(0L, nrow = 0, ncol = 2, dimnames = list(NULL, c("asm1", "asm2")))
  og_meta <- tibble::tibble(orthogroup_id = character(), size = integer())
  asm_meta <- tibble::tibble(assembly_name = c("asm1", "asm2"), n_orthogroups = c(0L, 0L))
  pa <- new_pa_matrix(mat, og_meta, asm_meta, type = "binary")

  expect_error(
    plot_upset(pa),
    regexp = "empty|orthogroup"
  )
})
