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

# Helper functions for score type tests ----------------------------------------

# Helper to create protein_collection with scores
make_scored_protein_collection <- function() {
  # asm1 has p1 (score=5), p2 (score=8), p3 (score=3)
  # asm2 has p1 (score=7)
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2", "asm1_p3"),
    sequence = c("AAA", "BBB", "CCC"),
    custom_score = c(5, 8, 3)
  )
  proteins_asm2 <- tibble::tibble(
    protein_id = c("asm2_p1"),
    sequence = c("DDD"),
    custom_score = c(7)
  )

  ps1 <- new_protein_set("asm1", proteins_asm1)
  ps2 <- new_protein_set("asm2", proteins_asm2)
  new_protein_collection(list(ps1, ps2))
}

# Orthogroup with paralogs having different scores
make_scored_orthogroup_result <- function() {
  # OG0001: asm1_p1 (score=5), asm1_p2 (score=8), asm2_p1 (score=7)
  # OG0002: asm1_p3 (score=3)
  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0001", "OG0002"),
    assembly = c("asm1", "asm1", "asm2", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2", "asm2_p1", "asm1_p3")
  )
  new_orthogroup_result(orthogroups, method = "test")
}

# Expected score matrix with max aggregation:
#        asm1  asm2
# OG0001    8     7   <- asm1 has max(5, 8)=8, asm2 has 7
# OG0002    3    NA   <- asm1 has 3, asm2 absent

# Expected score matrix with mean aggregation:
#        asm1  asm2
# OG0001  6.5     7   <- asm1 has mean(5, 8)=6.5
# OG0002    3    NA

# Expected score matrix with sum aggregation:
#        asm1  asm2
# OG0001   13     7   <- asm1 has sum(5, 8)=13
# OG0002    3    NA

# build_pa_matrix() score type tests -------------------------------------------

test_that("build_pa_matrix score type requires proteins argument", {
  ort <- make_scored_orthogroup_result()

  expect_error(
    build_pa_matrix(ort, type = "score", proteins = NULL),
    regexp = "proteins"
  )
})

test_that("build_pa_matrix score type looks up scores correctly", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score")

  # Single protein case: asm2_p1 has score 7

  expect_equal(pa$matrix["OG0001", "asm2"], 7)

  # Single protein case: asm1_p3 has score 3
  expect_equal(pa$matrix["OG0002", "asm1"], 3)
})

test_that("build_pa_matrix score aggregation max works", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        score_aggregation = "max")

  # OG0001 in asm1 has proteins with scores 5 and 8, max = 8
  expect_equal(pa$matrix["OG0001", "asm1"], 8)
})

test_that("build_pa_matrix score aggregation mean works", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        score_aggregation = "mean")

  # OG0001 in asm1 has proteins with scores 5 and 8, mean = 6.5
  expect_equal(pa$matrix["OG0001", "asm1"], 6.5)
})

test_that("build_pa_matrix score aggregation sum works", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        score_aggregation = "sum")

  # OG0001 in asm1 has proteins with scores 5 and 8, sum = 13
  expect_equal(pa$matrix["OG0001", "asm1"], 13)
})

test_that("build_pa_matrix score type uses NA for absent orthogroups", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score")

  # OG0002 is absent in asm2
  expect_true(is.na(pa$matrix["OG0002", "asm2"]))
})

test_that("build_pa_matrix score type sets type to 'score'", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score")

  expect_equal(pa$type, "score")
})

test_that("build_pa_matrix score type works with custom score column", {
  # Create protein collection with different score column name
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2"),
    sequence = c("AAA", "BBB"),
    my_score = c(10, 20)
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  proteins <- new_protein_collection(list(ps1))

  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001"),
    assembly = c("asm1", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        score_column = "my_score", score_aggregation = "max")

  expect_equal(pa$matrix["OG0001", "asm1"], 20)
})

# build_pa_matrix() score_threshold tests --------------------------------------

# Test data reminder (from make_scored_* helpers):
# OG0001: asm1_p1 (score=5), asm1_p2 (score=8), asm2_p1 (score=7) → max=8
# OG0002: asm1_p3 (score=3) → max=3
#
# With threshold=5:
# - OG0001 kept (max 8 >= 5)
# - OG0002 filtered out (max 3 < 5)

test_that("score_threshold requires proteins argument", {
  ort <- make_scored_orthogroup_result()

  expect_error(
    build_pa_matrix(ort, proteins = NULL, score_threshold = 5),
    regexp = "proteins"
  )
})

test_that("score_threshold filters orthogroups correctly", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, score_threshold = 5)

  # Only OG0001 should remain (max score 8 >= 5)
  # OG0002 filtered out (max score 3 < 5)
  expect_equal(nrow(pa$matrix), 1)
  expect_equal(rownames(pa$matrix), "OG0001")
})

test_that("score_threshold works with binary type", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "binary",
                        score_threshold = 5)

  # Binary matrix with only OG0001
  expect_equal(nrow(pa$matrix), 1)
  expect_equal(pa$matrix["OG0001", "asm1"], 1)
  expect_equal(pa$matrix["OG0001", "asm2"], 1)
  expect_equal(pa$type, "binary")
})

test_that("score_threshold works with count type", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "count",
                        score_threshold = 5)

  # Count matrix with only OG0001
  # OG0001 has 2 proteins in asm1, 1 in asm2
  expect_equal(nrow(pa$matrix), 1)
  expect_equal(pa$matrix["OG0001", "asm1"], 2)
  expect_equal(pa$matrix["OG0001", "asm2"], 1)
  expect_equal(pa$type, "count")
})

test_that("score_threshold works with score type", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, type = "score",
                        score_threshold = 5, score_aggregation = "max")

  # Score matrix with only OG0001
  expect_equal(nrow(pa$matrix), 1)
  expect_equal(pa$matrix["OG0001", "asm1"], 8)  # max(5, 8)
  expect_equal(pa$matrix["OG0001", "asm2"], 7)
  expect_equal(pa$type, "score")
})

test_that("score_threshold uses max score per orthogroup", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  # Threshold of 6: OG0001 max is 8 (passes), even though asm1_p1 has score 5
  pa <- build_pa_matrix(ort, proteins = proteins, score_threshold = 6)

  expect_equal(nrow(pa$matrix), 1)
  expect_equal(rownames(pa$matrix), "OG0001")
})

test_that("score_threshold stored in result", {
  ort <- make_scored_orthogroup_result()
  proteins <- make_scored_protein_collection()

  pa <- build_pa_matrix(ort, proteins = proteins, score_threshold = 5)

  expect_equal(pa$threshold, 5)
})

test_that("score_threshold with custom score_column", {
  # Create protein collection with different score column name
  proteins_asm1 <- tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2"),
    sequence = c("AAA", "BBB"),
    my_score = c(10, 2)  # p1 high, p2 low
  )
  ps1 <- new_protein_set("asm1", proteins_asm1)
  proteins <- new_protein_collection(list(ps1))

  orthogroups <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0002"),
    assembly = c("asm1", "asm1"),
    protein_id = c("asm1_p1", "asm1_p2")
  )
  ort <- new_orthogroup_result(orthogroups, method = "test")

  # Threshold 5 with my_score column
  pa <- build_pa_matrix(ort, proteins = proteins, score_threshold = 5,
                        score_column = "my_score")

  # Only OG0001 should remain (score 10 >= 5)
  expect_equal(nrow(pa$matrix), 1)
  expect_equal(rownames(pa$matrix), "OG0001")
})
