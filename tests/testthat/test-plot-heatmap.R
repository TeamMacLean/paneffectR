# Tests for plot_heatmap()

# Reuse helper from test-matrix.R
make_test_orthogroup_result <- function() {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002", "OG0002", "OG0003"),
    assembly = c("asm1", "asm2", "asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm2_p1", "asm1_p2", "asm2_p2", "asm1_p3")
  )
  new_orthogroup_result(orthogroups, method = "test")
}

# Helper for score type tests
make_scored_protein_collection <- function() {
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2", "asm1_p3"),
    sequence = c("AAA", "BBB", "CCC"),
    custom_score = c(5, 8, 3)
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = c("asm2_p1", "asm2_p2"),
    sequence = c("DDD", "EEE"),
    custom_score = c(7, 6)
  )

  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  new_protein_collection(list(ps1, ps2))
}

# plot_heatmap() basic tests ---------------------------------------------------

test_that("plot_heatmap returns ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())
  ht <- plot_heatmap(pa)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap validates input is pa_matrix", {
  skip_if_not_installed("ComplexHeatmap")

  expect_error(
    plot_heatmap("not a pa_matrix"),
    regexp = "pa_matrix"
  )
  expect_error(
    plot_heatmap(list(matrix = matrix(1:4, 2, 2))),
    regexp = "pa_matrix"
  )
  expect_error(
    plot_heatmap(NULL),
    regexp = "pa_matrix"
  )
})

test_that("plot_heatmap works with binary pa_matrix", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result(), type = "binary")
  ht <- plot_heatmap(pa)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap works with count pa_matrix", {
  skip_if_not_installed("ComplexHeatmap")

  # Create orthogroup with paralogs for count type
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001"),
    assembly = c("asm1", "asm1", "asm2"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, type = "count")

  ht <- plot_heatmap(pa)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap works with score pa_matrix", {
  skip_if_not_installed("ComplexHeatmap")

  ort <- make_test_orthogroup_result()
  proteins <- make_scored_protein_collection()
  pa <- build_pa_matrix(ort, proteins = proteins, type = "score")

  ht <- plot_heatmap(pa)

  expect_s4_class(ht, "Heatmap")
})

# plot_heatmap() clustering parameter tests ------------------------------------

test_that("plot_heatmap respects cluster_rows = TRUE", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())
  ht <- plot_heatmap(pa, cluster_rows = TRUE)

  expect_s4_class(ht, "Heatmap")
  # Verify clustering is enabled (slot check)
  expect_true(ht@row_dend_param$cluster)
})

test_that("plot_heatmap respects cluster_rows = FALSE", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())
  ht <- plot_heatmap(pa, cluster_rows = FALSE)

  expect_s4_class(ht, "Heatmap")
  # Verify clustering is disabled
  expect_false(ht@row_dend_param$cluster)
})

test_that("plot_heatmap respects cluster_cols = TRUE", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())
  ht <- plot_heatmap(pa, cluster_cols = TRUE)

  expect_s4_class(ht, "Heatmap")
  # Verify clustering is enabled
  expect_true(ht@column_dend_param$cluster)
})

test_that("plot_heatmap respects cluster_cols = FALSE", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())
  ht <- plot_heatmap(pa, cluster_cols = FALSE)

  expect_s4_class(ht, "Heatmap")
  # Verify clustering is disabled
  expect_false(ht@column_dend_param$cluster)
})

# plot_heatmap() edge cases ----------------------------------------------------

test_that("plot_heatmap handles single-row matrix", {
  skip_if_not_installed("ComplexHeatmap")

  # Single orthogroup
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001"),
    assembly = c("asm1", "asm2"),
    protein_id = c("asm1_p1", "asm2_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort)

  # Should not error, clustering automatically disabled for single row
  ht <- plot_heatmap(pa)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap handles single-column matrix", {
  skip_if_not_installed("ComplexHeatmap")

  # Single assembly
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    assembly = c("asm1", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort)

  # Should not error, clustering automatically disabled for single column
  ht <- plot_heatmap(pa)

  expect_s4_class(ht, "Heatmap")
})

# plot_heatmap() display parameter tests ---------------------------------------

test_that("plot_heatmap respects show_row_names parameter", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  ht_show <- plot_heatmap(pa, show_row_names = TRUE)
  ht_hide <- plot_heatmap(pa, show_row_names = FALSE)

  expect_true(ht_show@row_names_param$show)
  expect_false(ht_hide@row_names_param$show)
})

test_that("plot_heatmap respects show_col_names parameter", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  ht_show <- plot_heatmap(pa, show_col_names = TRUE)
  ht_hide <- plot_heatmap(pa, show_col_names = FALSE)

  expect_true(ht_show@column_names_param$show)
  expect_false(ht_hide@column_names_param$show)
})

test_that("plot_heatmap accepts custom color parameter", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())
  custom_colors <- c("gray90", "darkgreen")

  ht <- plot_heatmap(pa, color = custom_colors)

  expect_s4_class(ht, "Heatmap")
})
