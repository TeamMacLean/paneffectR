# Tests for build_pa_matrix() and related functions

# Helper function to create test orthogroup_result
make_test_orthogroup_result <- function() {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002", "OG0002", "OG0003"),
    assembly = c("asm1", "asm2", "asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm2_p1", "asm1_p2", "asm2_p2", "asm1_p3")
  )
  new_orthogroup_result(orthogroups, method = "test")
}

# Expected matrix for make_test_orthogroup_result():
#        asm1  asm2
# OG0001    1     1
# OG0002    1     1
# OG0003    1     0

# build_pa_matrix() tests -----------------------------------------------------

test_that("build_pa_matrix validates input is orthogroup_result", {
  expect_error(
    build_pa_matrix("not an orthogroup_result"),
    regexp = "orthogroup_result"
  )
  expect_error(
    build_pa_matrix(list(orthogroups = tibble::tibble())),
    regexp = "orthogroup_result"
  )
  expect_error(
    build_pa_matrix(NULL),
    regexp = "orthogroup_result"
  )
})

test_that("build_pa_matrix creates correct dimensions", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  # 3 orthogroups, 2 assemblies

  expect_equal(nrow(pa$matrix), 3)
  expect_equal(ncol(pa$matrix), 2)
})

test_that("build_pa_matrix produces binary values (0 or 1)",
 {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort, type = "binary")

  unique_values <- unique(as.vector(pa$matrix))
  expect_true(all(unique_values %in% c(0, 1)))
})

test_that("build_pa_matrix detects presence correctly", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  # OG0001 has proteins from both asm1 and asm2
  expect_equal(pa$matrix["OG0001", "asm1"], 1)
  expect_equal(pa$matrix["OG0001", "asm2"], 1)

  # OG0002 has proteins from both asm1 and asm2
  expect_equal(pa$matrix["OG0002", "asm1"], 1)
  expect_equal(pa$matrix["OG0002", "asm2"], 1)

  # OG0003 has protein only from asm1
  expect_equal(pa$matrix["OG0003", "asm1"], 1)
})

test_that("build_pa_matrix detects absence correctly", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  # OG0003 is absent in asm2
  expect_equal(pa$matrix["OG0003", "asm2"], 0)
})

test_that("build_pa_matrix row names match orthogroup IDs", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  expect_equal(sort(rownames(pa$matrix)), c("OG0001", "OG0002", "OG0003"))
})

test_that("build_pa_matrix column names match assembly names", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  expect_equal(sort(colnames(pa$matrix)), c("asm1", "asm2"))
})

test_that("build_pa_matrix returns pa_matrix class", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  expect_s3_class(pa, "pa_matrix")
})

test_that("build_pa_matrix orthogroups metadata has correct size", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  # Check orthogroups tibble has size column
  expect_true("size" %in% names(pa$orthogroups))

  # OG0001 has 2 proteins, OG0002 has 2 proteins, OG0003 has 1 protein
  og_sizes <- pa$orthogroups |>
    dplyr::arrange(.data$orthogroup_id)

  expect_equal(og_sizes$size[og_sizes$orthogroup_id == "OG0001"], 2)
  expect_equal(og_sizes$size[og_sizes$orthogroup_id == "OG0002"], 2)
  expect_equal(og_sizes$size[og_sizes$orthogroup_id == "OG0003"], 1)
})

test_that("build_pa_matrix assemblies metadata has correct n_orthogroups", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  # Check assemblies tibble has n_orthogroups column
  expect_true("n_orthogroups" %in% names(pa$assemblies))

  # asm1 is in all 3 orthogroups, asm2 is in 2 orthogroups
  asm_data <- pa$assemblies

  expect_equal(asm_data$n_orthogroups[asm_data$assembly_name == "asm1"], 3)
  expect_equal(asm_data$n_orthogroups[asm_data$assembly_name == "asm2"], 2)
})

test_that("build_pa_matrix sets type to 'binary' by default", {
  ort <- make_test_orthogroup_result()
  pa <- build_pa_matrix(ort)

  expect_equal(pa$type, "binary")
})

test_that("build_pa_matrix handles orthogroups with paralogs", {
  # Create orthogroup with multiple proteins from same assembly (paralogs)
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001"),
    assembly = c("asm1", "asm1", "asm2"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  pa <- build_pa_matrix(ort, type = "binary")

  # Even with paralogs, binary should still be 1 (not 2)
  expect_equal(pa$matrix["OG0001", "asm1"], 1)
  expect_equal(pa$matrix["OG0001", "asm2"], 1)
})

# Helper function for count type tests with paralogs
make_paralog_orthogroup_result <- function() {
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001", "OG0002"),
    assembly = c("asm1", "asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1", "asm1_p3")
  )
  new_orthogroup_result(orthogroups, method = "test")
}
# Expected count matrix for make_paralog_orthogroup_result():
#        asm1  asm2
# OG0001    2     1   <- asm1 has 2 paralogs
# OG0002    1     0

# build_pa_matrix() count type tests ---------------------------------------------

test_that("build_pa_matrix count type produces integer counts", {
  ort <- make_paralog_orthogroup_result()
  pa <- build_pa_matrix(ort, type = "count")

  # All values should be non-negative integers
  expect_true(all(pa$matrix >= 0))
  expect_true(all(pa$matrix == as.integer(pa$matrix)))
})

test_that("build_pa_matrix count type counts paralogs correctly", {
  ort <- make_paralog_orthogroup_result()
  pa <- build_pa_matrix(ort, type = "count")

  # OG0001 has 2 proteins from asm1 (paralogs)
  expect_equal(pa$matrix["OG0001", "asm1"], 2)
})

test_that("build_pa_matrix count type single copy remains 1", {
  ort <- make_paralog_orthogroup_result()
  pa <- build_pa_matrix(ort, type = "count")

  # OG0001 has 1 protein from asm2
  expect_equal(pa$matrix["OG0001", "asm2"], 1)

  # OG0002 has 1 protein from asm1
  expect_equal(pa$matrix["OG0002", "asm1"], 1)
})

test_that("build_pa_matrix count type absence is 0", {
  ort <- make_paralog_orthogroup_result()
  pa <- build_pa_matrix(ort, type = "count")

  # OG0002 has no proteins from asm2
  expect_equal(pa$matrix["OG0002", "asm2"], 0)
})

test_that("build_pa_matrix count type sets type to 'count'", {
  ort <- make_paralog_orthogroup_result()
  pa <- build_pa_matrix(ort, type = "count")

  expect_equal(pa$type, "count")
})
