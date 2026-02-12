# Tests for singleton utility functions

# Helper function to create test orthogroup_result with singletons
make_orthogroup_with_singletons <- function() {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002", "OG0002"),
    assembly = c("asm1", "asm2", "asm1", "asm2"),
    protein_id = c("asm1_p1", "asm2_p1", "asm1_p2", "asm2_p2")
  )
  singletons <- tibble::tibble(
    assembly = c("asm1", "asm2", "asm2"),
    protein_id = c("asm1_single1", "asm2_single1", "asm2_single2")
  )
  new_orthogroup_result(orthogroups, method = "test", singletons = singletons)
}

# get_singletons() tests -----------------------------------------------------

test_that("get_singletons returns singletons tibble", {
  ort <- make_orthogroup_with_singletons()

  result <- get_singletons(ort)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
})

test_that("get_singletons includes protein_id and assembly columns", {
  ort <- make_orthogroup_with_singletons()

  result <- get_singletons(ort)

  expect_true("protein_id" %in% names(result))
  expect_true("assembly" %in% names(result))
})

test_that("get_singletons returns correct protein_ids", {
  ort <- make_orthogroup_with_singletons()

  result <- get_singletons(ort)

  expect_equal(sort(result$protein_id), c("asm1_single1", "asm2_single1", "asm2_single2"))
})

test_that("get_singletons returns empty tibble when no singletons", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001"),
    assembly = c("asm1"),
    protein_id = c("asm1_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  result <- get_singletons(ort)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_true("protein_id" %in% names(result))
  expect_true("assembly" %in% names(result))
})

test_that("get_singletons validates input is orthogroup_result", {
  expect_error(
    get_singletons("not an orthogroup_result"),
    regexp = "orthogroup_result"
  )
})

# n_singletons() tests -------------------------------------------------------

test_that("n_singletons returns correct count", {
  ort <- make_orthogroup_with_singletons()

  result <- n_singletons(ort)

  expect_equal(result, 3)
})

test_that("n_singletons returns 0 when no singletons", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001"),
    assembly = c("asm1"),
    protein_id = c("asm1_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  result <- n_singletons(ort)

  expect_equal(result, 0)
})

test_that("n_singletons validates input is orthogroup_result", {
  expect_error(
    n_singletons("not an orthogroup_result"),
    regexp = "orthogroup_result"
  )
})

# singletons_by_assembly() tests ---------------------------------------------

test_that("singletons_by_assembly returns tibble with expected columns", {
  ort <- make_orthogroup_with_singletons()

  result <- singletons_by_assembly(ort)

  expect_s3_class(result, "tbl_df")
  expect_true("assembly" %in% names(result))
  expect_true("n_singletons" %in% names(result))
})

test_that("singletons_by_assembly counts correctly per assembly", {
  ort <- make_orthogroup_with_singletons()

  result <- singletons_by_assembly(ort)

  # asm1 has 1 singleton, asm2 has 2 singletons
  asm1_count <- result$n_singletons[result$assembly == "asm1"]
  asm2_count <- result$n_singletons[result$assembly == "asm2"]

  expect_equal(asm1_count, 1)
  expect_equal(asm2_count, 2)
})

test_that("singletons_by_assembly calculates percentage when proteins provided", {
  ort <- make_orthogroup_with_singletons()

  # Create protein_collection with known protein counts
  # asm1: 3 proteins total (p1, p2, single1) -> 1 singleton -> 33.33%
  # asm2: 4 proteins total (p1, p2, single1, single2) -> 2 singletons -> 50%
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2", "asm1_single1"),
    sequence = c("AAA", "BBB", "CCC")
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = c("asm2_p1", "asm2_p2", "asm2_single1", "asm2_single2"),
    sequence = c("DDD", "EEE", "FFF", "GGG")
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  proteins <- new_protein_collection(list(ps1, ps2))

  result <- singletons_by_assembly(ort, proteins = proteins)

  expect_true("n_total" %in% names(result))
  expect_true("pct_singleton" %in% names(result))

  # Check percentages - use [[1]] to extract scalar values
  asm1_row <- result[result$assembly == "asm1", ]
  asm2_row <- result[result$assembly == "asm2", ]

  expect_equal(asm1_row$n_total[[1]], 3)
  expect_equal(asm1_row$pct_singleton[[1]], 1/3 * 100, tolerance = 0.01)

  expect_equal(asm2_row$n_total[[1]], 4)
  expect_equal(asm2_row$pct_singleton[[1]], 50)
})

test_that("singletons_by_assembly works without proteins argument", {
  ort <- make_orthogroup_with_singletons()

  result <- singletons_by_assembly(ort)

  # Should work but without n_total and pct_singleton columns
  expect_true("assembly" %in% names(result))
  expect_true("n_singletons" %in% names(result))
})

test_that("singletons_by_assembly handles empty singletons", {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001"),
    assembly = c("asm1"),
    protein_id = c("asm1_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  result <- singletons_by_assembly(ort)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("singletons_by_assembly validates input is orthogroup_result", {
  expect_error(
    singletons_by_assembly("not an orthogroup_result"),
    regexp = "orthogroup_result"
  )
})
