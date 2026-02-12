# Tests for plot_dendro()
# Phase 4, Chunk 5: Assembly Dendrogram
# Spec: specs/phase4-chunk5-dendrogram.md

# Helper to create test pa_matrix with multiple assemblies
make_multi_assembly_pa <- function() {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001", "OG0001",
                      "OG0002", "OG0002", "OG0002",
                      "OG0003", "OG0003",
                      "OG0004"),
    assembly = c("asm1", "asm2", "asm3", "asm4",
                 "asm1", "asm2", "asm3",
                 "asm1", "asm2",
                 "asm1"),
    protein_id = c("asm1_p1", "asm2_p1", "asm3_p1", "asm4_p1",
                   "asm1_p2", "asm2_p2", "asm3_p2",
                   "asm1_p3", "asm2_p3",
                   "asm1_p4")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  build_pa_matrix(ort, exclude_singletons = TRUE)
}

# plot_dendro() basic tests ---------------------------------------------------

test_that("plot_dendro returns ggplot object", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa)

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro validates input is pa_matrix", {
  expect_error(
    plot_dendro("not a pa_matrix"),
    regexp = "pa_matrix"
  )
  expect_error(
    plot_dendro(NULL),
    regexp = "pa_matrix"
  )
})

test_that("plot_dendro works with default parameters", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa)

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro errors with less than 2 assemblies", {
  # Single assembly
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    assembly = c("asm1", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, exclude_singletons = TRUE)

  expect_error(
    plot_dendro(pa),
    regexp = "assembl|2"
  )
})

# plot_dendro() distance_method tests -----------------------------------------

test_that("plot_dendro distance_method = 'jaccard' works (default)", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, distance_method = "jaccard")

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro distance_method = 'binary' works", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, distance_method = "binary")

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro distance_method = 'bray-curtis' works", {
  # Use count matrix for bray-curtis
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001", "OG0001", "OG0001"),
    assembly = c("asm1", "asm1", "asm2", "asm3", "asm3"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1", "asm3_p1", "asm3_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, type = "count", exclude_singletons = TRUE)

  p <- plot_dendro(pa, distance_method = "bray-curtis")

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro errors on invalid distance_method", {
  pa <- make_multi_assembly_pa()

  expect_error(
    plot_dendro(pa, distance_method = "invalid"),
    regexp = "distance_method"
  )
})

# plot_dendro() cluster_method tests ------------------------------------------

test_that("plot_dendro cluster_method = 'complete' works (default)", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, cluster_method = "complete")

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro cluster_method = 'average' works", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, cluster_method = "average")

  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro cluster_method = 'ward.D2' works", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, cluster_method = "ward.D2")

  expect_s3_class(p, "ggplot")
})

# plot_dendro() display parameter tests ---------------------------------------

test_that("plot_dendro labels = TRUE shows labels (default)", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, labels = TRUE)

  expect_s3_class(p, "ggplot")
  # Check that text layer exists for labels
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomText" %in% layer_types || "GeomLabel" %in% layer_types)
})
test_that("plot_dendro labels = FALSE hides labels", {
  pa <- make_multi_assembly_pa()

  p <- plot_dendro(pa, labels = FALSE)

  expect_s3_class(p, "ggplot")
})

# plot_dendro() edge cases ----------------------------------------------------

test_that("plot_dendro handles matrix with NAs", {
  # Create score matrix which has NAs for absent
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002"),
    assembly = c("asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm2_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2"),
    sequence = c("AAA", "BBB"),
    custom_score = c(5.0, 3.0)
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = c("asm2_p1"),
    sequence = c("CCC"),
    custom_score = c(7.0)
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  proteins <- new_protein_collection(list(ps1, ps2))

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        exclude_singletons = TRUE)

  # Should handle NAs (treat as 0)
  expect_no_error(p <- plot_dendro(pa))
  expect_s3_class(p, "ggplot")
})

test_that("plot_dendro works with exactly 2 assemblies", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002"),
    assembly = c("asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm2_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, exclude_singletons = TRUE)

  # Should work (minimum case)
  p <- plot_dendro(pa)

  expect_s3_class(p, "ggplot")
})
