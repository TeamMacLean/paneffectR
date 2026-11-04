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

# =============================================================================
# Chunk 2: Annotation and Distance Metric Tests
# =============================================================================

# plot_heatmap() row_annotation tests -----------------------------------------

test_that("plot_heatmap works with row_annotation (single column)", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  row_ann <- data.frame(
    size = c(2, 2, 1),
    row.names = c("OG0001", "OG0002", "OG0003")
  )

  ht <- plot_heatmap(pa, row_annotation = row_ann)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap works with row_annotation (multiple columns)", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  row_ann <- data.frame(
    size = c(2, 2, 1),
    category = c("shared", "shared", "unique"),
    row.names = c("OG0001", "OG0002", "OG0003")
  )

  ht <- plot_heatmap(pa, row_annotation = row_ann)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap works with partial row_annotation match", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  # Only annotate 2 of 3 orthogroups
  row_ann <- data.frame(
    size = c(2, 2),
    row.names = c("OG0001", "OG0002")
  )

  # Should work with warning about partial match
  expect_warning(
    ht <- plot_heatmap(pa, row_annotation = row_ann),
    regexp = "annotation|match"
  )

  expect_s4_class(ht, "Heatmap")
})

# plot_heatmap() col_annotation tests -----------------------------------------

test_that("plot_heatmap works with col_annotation", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  col_ann <- data.frame(
    species = c("P. destructor", "P. destructor"),
    row.names = c("asm1", "asm2")
  )

  ht <- plot_heatmap(pa, col_annotation = col_ann)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap works with both row and col annotations", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  row_ann <- data.frame(
    size = c(2, 2, 1),
    row.names = c("OG0001", "OG0002", "OG0003")
  )

  col_ann <- data.frame(
    species = c("P. destructor", "P. destructor"),
    row.names = c("asm1", "asm2")
  )

  ht <- plot_heatmap(pa, row_annotation = row_ann, col_annotation = col_ann)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap warns on mismatched col_annotation names", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  # Only annotate 1 of 2 assemblies
  col_ann <- data.frame(
    species = c("P. destructor"),
    row.names = c("asm1")
  )

  expect_warning(
    ht <- plot_heatmap(pa, col_annotation = col_ann),
    regexp = "annotation|match"
  )

  expect_s4_class(ht, "Heatmap")
})

# plot_heatmap() distance_method tests ----------------------------------------

test_that("plot_heatmap distance_method = 'binary' works (default)", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  ht <- plot_heatmap(pa, cluster_rows = TRUE, distance_method = "binary")

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap distance_method = 'jaccard' works", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  ht <- plot_heatmap(pa, cluster_rows = TRUE, distance_method = "jaccard")

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap distance_method = 'bray-curtis' works", {
  skip_if_not_installed("ComplexHeatmap")

  # Use count matrix for bray-curtis
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001", "OG0002"),
    assembly = c("asm1", "asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1", "asm1_p3")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")
  pa <- build_pa_matrix(ort, type = "count")

  ht <- plot_heatmap(pa, cluster_rows = TRUE, distance_method = "bray-curtis")

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_heatmap errors on invalid distance_method", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- build_pa_matrix(make_test_orthogroup_result())

  expect_error(
    plot_heatmap(pa, distance_method = "invalid_method"),
    regexp = "distance_method"
  )
})

# plot_heatmap styles ----------------------------------------------------------

visual_pa <- function(name) {
  readRDS(file.path(system.file("testdata", "visual", package = "paneffectR"), name))
}

test_that("plot_heatmap accepts the three styles for every matrix type", {
  skip_if_not_installed("ComplexHeatmap")

  pa_bin <- visual_pa("pa_binary.rds")
  pa_score <- visual_pa("pa_score.rds")
  pa_count <- build_pa_matrix(visual_pa("clusters_visual.rds"), type = "count")

  for (s in c("clean", "pan", "mono")) {
    for (pa in list(pa_bin, pa_score, pa_count)) {
      ht <- plot_heatmap(pa, style = s)
      expect_s4_class(ht, "Heatmap")
      # Drawing exercises ComplexHeatmap's own argument checks
      expect_no_error(withr::with_pdf(NULL, ComplexHeatmap::draw(ht)))
    }
  }
})

test_that("plot_heatmap rejects an unknown style", {
  skip_if_not_installed("ComplexHeatmap")
  expect_error(plot_heatmap(visual_pa("pa_binary.rds"), style = "fancy"), "style")
})

test_that("style legends are labelled Absent/Present with a meaningful title", {
  skip_if_not_installed("ComplexHeatmap")

  ht <- plot_heatmap(visual_pa("pa_binary.rds"))
  expect_equal(ht@name, "Orthogroup")
  expect_equal(ht@matrix_legend_param$labels, c("Absent", "Present"))

  ht <- plot_heatmap(visual_pa("pa_score.rds"))
  expect_equal(ht@name, "Score")
})

test_that("count matrices with a maximum of 1 use the binary colours", {
  skip_if_not_installed("ComplexHeatmap")

  pa_count <- build_pa_matrix(visual_pa("clusters_visual.rds"), type = "count")
  expect_equal(max(pa_count$matrix), 1)
  ht <- plot_heatmap(pa_count)
  expect_equal(ht@matrix_legend_param$labels, c("Absent", "Present"))
})

test_that("classify_rows labels Core, Accessory, Unique and Absent", {
  m <- rbind(c(1, 1, 1), c(1, 0, 1), c(0, 0, 1), c(0, 0, 0))
  cls <- classify_rows(m)
  expect_equal(as.character(cls), c("Core", "Accessory", "Unique", "Absent"))
  expect_equal(levels(cls), c("Core", "Accessory", "Unique", "Absent"))

  # unused levels are dropped
  expect_equal(levels(classify_rows(m[1:2, , drop = FALSE])), c("Core", "Accessory"))

  # a single column: every present row is Core
  expect_equal(as.character(classify_rows(matrix(c(1, 1, 0), ncol = 1))),
               c("Core", "Core", "Absent"))
})

test_that("pan style splits rows by class and colours binary cells by class", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- visual_pa("pa_binary.rds")
  ht <- ComplexHeatmap::draw(plot_heatmap(pa, style = "pan"))
  ro <- ComplexHeatmap::row_order(ht)
  expect_setequal(names(ro), c("Core", "Accessory", "Unique"))
  expect_equal(sum(lengths(ro)), nrow(pa$matrix))

  # every row in the Core slice is present in all assemblies
  expect_true(all(rowSums(pa$matrix[ro$Core, , drop = FALSE]) == ncol(pa$matrix)))
  expect_true(all(rowSums(pa$matrix[ro$Unique, , drop = FALSE]) == 1))

  # cells are coded by class; legend lists the classes
  lp <- plot_heatmap(pa, style = "pan")@matrix_legend_param
  expect_equal(lp$labels, c("Core", "Accessory", "Unique", "Absent"))
})

test_that("pan style keeps the pattern-ordering fallback for large matrices", {
  skip_if_not_installed("ComplexHeatmap")

  set.seed(2)
  n <- 300
  m <- matrix(rbinom(n * 4, 1, 0.5), nrow = n,
              dimnames = list(sprintf("OG%04d", seq_len(n)), paste0("asm", 1:4)))
  storage.mode(m) <- "integer"
  pa <- new_pa_matrix(
    matrix = m,
    orthogroups = tibble::tibble(orthogroup_id = rownames(m), size = rowSums(m)),
    assemblies = tibble::tibble(assembly_name = colnames(m), n_orthogroups = colSums(m))
  )

  expect_message(
    ht <- plot_heatmap(pa, style = "pan", max_cluster_rows = 100),
    "distinct presence/absence pattern"
  )
  ht <- withr::with_pdf(NULL, ComplexHeatmap::draw(ht))
  ro <- ComplexHeatmap::row_order(ht)
  # identical patterns are adjacent within the Accessory slice
  acc <- ro$Accessory
  keys <- apply(m[acc, , drop = FALSE], 1, paste, collapse = "")
  expect_equal(length(rle(keys)$values), length(unique(keys)))
})

test_that("pan style works with row and column annotations", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- visual_pa("pa_binary.rds")
  row_ann <- data.frame(row.names = rownames(pa$matrix),
                        grp = rep(c("x", "y"), length.out = nrow(pa$matrix)))
  col_ann <- data.frame(row.names = colnames(pa$matrix),
                        host = rep(c("a", "b"), length.out = ncol(pa$matrix)))
  ht <- plot_heatmap(pa, style = "pan", row_annotation = row_ann, col_annotation = col_ann)
  expect_no_error(withr::with_pdf(NULL, ComplexHeatmap::draw(ht)))
})

test_that("user color overrides the style palette", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- visual_pa("pa_binary.rds")
  for (s in c("clean", "pan", "mono")) {
    ht <- plot_heatmap(pa, style = s, color = c("grey95", "purple"))
    expect_no_error(withr::with_pdf(NULL, ComplexHeatmap::draw(ht)))
  }
})

test_that("arguments in ... override style settings", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- visual_pa("pa_binary.rds")
  ht <- plot_heatmap(pa, column_names_rot = 90, name = "Custom",
                     rect_gp = grid::gpar(col = NA))
  expect_equal(ht@name, "Custom")
  expect_equal(ht@column_names_param$rot, 90)
  expect_no_error(withr::with_pdf(NULL, ComplexHeatmap::draw(ht)))
})

test_that("styles handle single-row, single-column and constant-score matrices", {
  skip_if_not_installed("ComplexHeatmap")

  pa <- visual_pa("pa_binary.rds")
  ps <- visual_pa("pa_score.rds")
  ps$matrix[!is.na(ps$matrix)] <- 5

  for (s in c("clean", "pan", "mono")) {
    expect_s4_class(plot_heatmap(pa[1, ], style = s), "Heatmap")
    expect_s4_class(plot_heatmap(pa[, 1], style = s), "Heatmap")
    expect_s4_class(plot_heatmap(ps, style = s), "Heatmap")
  }
})

test_that("grid lines are dropped above 200 rows", {
  skip_if_not_installed("ComplexHeatmap")

  set.seed(3)
  n <- 250
  m <- matrix(rbinom(n * 3, 1, 0.5), nrow = n,
              dimnames = list(sprintf("OG%04d", seq_len(n)), paste0("asm", 1:3)))
  storage.mode(m) <- "integer"
  pa <- new_pa_matrix(
    matrix = m,
    orthogroups = tibble::tibble(orthogroup_id = rownames(m), size = rowSums(m)),
    assemblies = tibble::tibble(assembly_name = colnames(m), n_orthogroups = colSums(m))
  )
  ht <- plot_heatmap(pa)
  expect_true(is.na(ht@matrix_param$gp$col))
  ht_small <- plot_heatmap(pa[1:50, ])
  expect_equal(ht_small@matrix_param$gp$col, "white")
})
