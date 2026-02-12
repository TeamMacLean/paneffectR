# Tests for plot_scores()
# Phase 4, Chunk 4: Score Distributions
# Spec: specs/phase4-chunk4-score-distributions.md

# Helper to create test protein_collection with scores
make_scored_proteins <- function() {
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2", "asm1_p3"),
    sequence = c("AAA", "BBB", "CCC"),
    custom_score = c(5.0, 8.0, 3.0)
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = c("asm2_p1", "asm2_p2"),
    sequence = c("DDD", "EEE"),
    custom_score = c(7.0, 6.0)
  )

  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  new_protein_collection(list(ps1, ps2))
}

# plot_scores() basic tests ---------------------------------------------------

test_that("plot_scores returns ggplot object", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins)

  expect_s3_class(p, "ggplot")
})

test_that("plot_scores validates input is protein_collection", {
  expect_error(
    plot_scores("not a protein_collection"),
    regexp = "protein_collection"
  )
  expect_error(
    plot_scores(NULL),
    regexp = "protein_collection"
  )
})

test_that("plot_scores works with default parameters", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins)

  expect_s3_class(p, "ggplot")
})

test_that("plot_scores errors when score_column doesn't exist", {
  proteins <- make_scored_proteins()

  expect_error(
    plot_scores(proteins, score_column = "nonexistent_column"),
    regexp = "score|column|found"
  )
})

# plot_scores() parameter tests -----------------------------------------------

test_that("plot_scores by_assembly = TRUE creates faceted plot", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins, by_assembly = TRUE)

  expect_s3_class(p, "ggplot")
  # Check that faceting is present
  expect_true("FacetWrap" %in% class(p$facet) || "FacetGrid" %in% class(p$facet))
})

test_that("plot_scores by_assembly = FALSE creates single plot", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins, by_assembly = FALSE)

  expect_s3_class(p, "ggplot")
  # Check that faceting is NOT present (FacetNull)
  expect_true("FacetNull" %in% class(p$facet))
})

test_that("plot_scores threshold adds vertical line", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins, threshold = 5.0)

  expect_s3_class(p, "ggplot")
  # Check that vline layer exists
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomVline" %in% layer_types)
})

test_that("plot_scores plot_type = 'density' creates density plot", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins, plot_type = "density")

  expect_s3_class(p, "ggplot")
  # Check for density geom
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomDensity" %in% layer_types || "GeomArea" %in% layer_types)
})

test_that("plot_scores plot_type = 'histogram' creates histogram", {
  proteins <- make_scored_proteins()

  p <- plot_scores(proteins, plot_type = "histogram")

  expect_s3_class(p, "ggplot")
  # Check for histogram/bar geom
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomBar" %in% layer_types || "GeomHistogram" %in% layer_types)
})

test_that("plot_scores errors on invalid plot_type", {
  proteins <- make_scored_proteins()

  expect_error(
    plot_scores(proteins, plot_type = "invalid"),
    regexp = "plot_type"
  )
})

test_that("plot_scores custom score_column works", {
  # Create proteins with different score column
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2"),
    sequence = c("AAA", "BBB"),
    my_score = c(10.0, 20.0)
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  proteins <- new_protein_collection(list(ps1))

  p <- plot_scores(proteins, score_column = "my_score")

  expect_s3_class(p, "ggplot")
})

# plot_scores() edge cases ----------------------------------------------------

test_that("plot_scores handles NAs in score column", {
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2", "asm1_p3"),
    sequence = c("AAA", "BBB", "CCC"),
    custom_score = c(5.0, NA, 3.0)
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  proteins <- new_protein_collection(list(ps1))

  # Should handle NAs gracefully (remove them)
  expect_no_error(p <- plot_scores(proteins))
  expect_s3_class(p, "ggplot")
})

test_that("plot_scores works when only some assemblies have scores", {
  # asm1 has scores, asm2 doesn't
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1"),
    sequence = c("AAA"),
    custom_score = c(5.0)
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = c("asm2_p1"),
    sequence = c("BBB")
    # No custom_score column
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  proteins <- new_protein_collection(list(ps1, ps2))

  # Should work with just the assemblies that have scores
  expect_no_error(p <- plot_scores(proteins))
  expect_s3_class(p, "ggplot")
})

test_that("plot_scores errors when no assemblies have scores", {
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1"),
    sequence = c("AAA")
    # No score column
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  proteins <- new_protein_collection(list(ps1))

  expect_error(
    plot_scores(proteins),
    regexp = "score|found"
  )
})
