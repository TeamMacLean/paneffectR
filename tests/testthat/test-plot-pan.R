# Tests for plot-pan.R - pan-genome structure visualizations

test_that("plot_pan_structure returns ggplot object", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  p <- plot_pan_structure(clusters)

  expect_s3_class(p, "ggplot")
})
test_that("plot_pan_structure works with exclude_singletons", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  p_with <- plot_pan_structure(clusters, exclude_singletons = FALSE)
  p_without <- plot_pan_structure(clusters, exclude_singletons = TRUE)

  expect_s3_class(p_with, "ggplot")
  expect_s3_class(p_without, "ggplot")

  # Check that singleton category is absent when excluded
  data_with <- ggplot2::ggplot_build(p_with)$data[[1]]
  data_without <- ggplot2::ggplot_build(p_without)$data[[1]]

  expect_gt(nrow(data_with), nrow(data_without))
})

test_that("plot_pan_structure accepts custom thresholds", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  # Should not error with custom thresholds
  p <- plot_pan_structure(
    clusters,
    rare_threshold = 0.20,
    accessory_threshold = 0.70
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_pan_structure accepts custom colors", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  custom_colors <- c(
    "Core" = "red",
    "Accessory" = "blue",
    "Rare" = "green",
    "Unique" = "yellow",
    "Singleton" = "purple"
  )

  p <- plot_pan_structure(clusters, colors = custom_colors)

  expect_s3_class(p, "ggplot")
})

test_that("plot_pan_structure validates input", {
  expect_error(
    plot_pan_structure("not a cluster"),
    "must be an.*orthogroup_result"
  )
})

test_that("plot_assembly_composition returns ggplot object", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  p <- plot_assembly_composition(clusters)

  expect_s3_class(p, "ggplot")
})

test_that("plot_assembly_composition works with position fill", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  p_stack <- plot_assembly_composition(clusters, position = "stack")
  p_fill <- plot_assembly_composition(clusters, position = "fill")

  expect_s3_class(p_stack, "ggplot")
  expect_s3_class(p_fill, "ggplot")
})

test_that("plot_assembly_composition works with exclude_singletons", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  p <- plot_assembly_composition(clusters, exclude_singletons = TRUE)

  expect_s3_class(p, "ggplot")
})

test_that("plot_assembly_composition validates input", {
  expect_error(
    plot_assembly_composition(list(not = "valid")),
    "must be an.*orthogroup_result"
  )
})

test_that("categorize_orthogroups returns expected structure", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))

  categories <- paneffectR:::categorize_orthogroups(clusters)

  expect_s3_class(categories, "tbl_df")
  expect_true(all(c("orthogroup_id", "n_present", "n_assemblies", "category") %in% names(categories)))
  expect_true(all(categories$category %in% c("Core", "Accessory", "Rare", "Unique")))
})
