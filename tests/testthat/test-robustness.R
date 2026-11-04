# Tests for input robustness and UX fixes from the new-user test pass (ISSUES.md)

get_testdata_path <- function(file = NULL) {
  path <- system.file("testdata", package = "paneffectR")
  if (path == "") {
    path <- test_path("../../inst/testdata")
  }
  if (!is.null(file)) {
    path <- file.path(path, file)
  }
  path
}

# Copy the three bundled assemblies (and companions) into a temp dir
copy_testdata <- function(dir = withr::local_tempdir(.local_envir = parent.frame())) {
  files <- list.files(get_testdata_path(), pattern = "^assembly[0-9]", full.names = TRUE)
  file.copy(files, dir)
  dir
}

# Add a free-text description to every header in a FASTA file
add_header_descriptions <- function(path, desc = "hypothetical protein [Pyricularia oryzae]") {
  lines <- readLines(path)
  is_header <- startsWith(lines, ">")
  lines[is_header] <- paste(lines[is_header], desc)
  writeLines(lines, path)
  invisible(path)
}

get_diamond_path <- function() {
  env_path <- Sys.getenv("PANEFFECTR_DIAMOND", unset = "")
  if (nzchar(env_path) && file.exists(env_path)) return(env_path)
  for (rel in c("../../this_project_env/bin/diamond",
                "../../../this_project_env/bin/diamond")) {
    local_path <- test_path(rel)
    if (file.exists(local_path)) return(normalizePath(local_path))
  }
  sys_path <- Sys.which("diamond")
  if (nzchar(sys_path)) return(unname(sys_path))
  NULL
}

# FASTA header handling (issue 2) ----------------------------------------------

test_that("load_fasta uses the first header token as protein_id and keeps the description", {
  tf <- withr::local_tempfile(fileext = ".faa")
  writeLines(c(
    ">XP_003708685.1 cytochrome P450 3A19 [Pyricularia oryzae 70-15]",
    "MKVLAAGIVALLLAAGCSSA",
    ">XP_000000002.1",
    "MSTNPKPQRKTKRNTNRRPQ",
    ">sp|P12345|TEST_HUMAN\tTab separated description",
    "MAAAAAAAAAGGGGG"
  ), tf)

  result <- load_fasta(tf)

  expect_named(result, c("protein_id", "sequence", "description"))
  expect_equal(result$protein_id, c("XP_003708685.1", "XP_000000002.1", "sp|P12345|TEST_HUMAN"))
  expect_equal(result$description,
               c("cytochrome P450 3A19 [Pyricularia oryzae 70-15]", NA, "Tab separated description"))
  expect_false(any(grepl("\\s", result$protein_id)))
})

test_that("load_fasta strips trailing stop codons and upper-cases sequences", {
  tf <- withr::local_tempfile(fileext = ".faa")
  writeLines(c(">p1", "mkvlaagiv*", ">p2", "MSTNPKPQR**", ">p3", "MAAAGG"), tf)

  result <- load_fasta(tf)

  expect_equal(result$sequence, c("MKVLAAGIV", "MSTNPKPQR", "MAAAGG"))
})

test_that("new_protein_set rejects protein IDs containing whitespace", {
  expect_error(
    new_protein_set("asm", tibble::tibble(protein_id = "p1 description", sequence = "MKV")),
    "whitespace"
  )
})

test_that("duplicate ID error names the assembly and the IDs", {
  expect_error(
    new_protein_set("asmX", tibble::tibble(protein_id = c("p1", "p1"), sequence = c("MKV", "MKV"))),
    "asmX.*p1"
  )
})

test_that("headers with descriptions cluster correctly end to end", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")

  d <- copy_testdata()
  for (f in list.files(d, pattern = "\\.faa$", full.names = TRUE)) add_header_descriptions(f)

  p <- load_proteins(d, score_dir = d, quiet = TRUE)
  expect_true(all(!is.na(unlist(lapply(p$assemblies, function(ps) ps$proteins$description)))))

  cl <- cluster_proteins(p, tool_path = diamond_path, quiet = TRUE)
  expect_false(anyNA(cl$orthogroups$assembly))
  expect_lt(n_singletons(cl), sum(p$summary$n_proteins))
  expect_gt(nrow(cl$orthogroups), 0)

  pa <- build_pa_matrix(cl, proteins = p, score_threshold = 5)
  expect_s3_class(pa, "pa_matrix")
  expect_equal(sort(colnames(pa$matrix)), c("assembly1", "assembly2", "assembly3"))
})

# Nucleotide input (issue 3) ---------------------------------------------------

test_that("nucleotide_fraction distinguishes protein from nucleotide", {
  expect_gt(nucleotide_fraction(c("ACGTACGTACGTNNACGT", "acgtacgt")), 0.95)
  expect_lt(nucleotide_fraction(c("MKVLAAGIVALLLAAGCSSA", "MSTNPKPQRKTKRNTNRRPQ")), 0.5)
  expect_true(is.na(nucleotide_fraction(character())))
  expect_true(is.na(nucleotide_fraction("")))
})

test_that("load_fasta rejects nucleotide FASTA with a clear message", {
  tf <- withr::local_tempfile(fileext = ".faa")
  writeLines(c(">gene1", "ATGGCTAGCTAGCTAGCTAGCTAGGATCGATCGATCG", ">gene2", "ATGCCCGGGTTTAAACCCGGG"), tf)
  expect_error(load_fasta(tf), "nucleotide")
})

test_that("cluster_proteins rejects nucleotide-like collections before calling DIAMOND", {
  ps1 <- new_protein_set("a", tibble::tibble(protein_id = "a1", sequence = "ATGCATGCATGCATGC"))
  ps2 <- new_protein_set("b", tibble::tibble(protein_id = "b1", sequence = "ATGCATGCATGCATGA"))
  pc <- new_protein_collection(list(ps1, ps2))
  # Error must come from validation, not from tool lookup
  expect_error(cluster_proteins(pc, tool_path = "/nonexistent/diamond"), "nucleotide")
})

test_that("cluster_proteins validates timeout", {
  ps1 <- new_protein_set("a", tibble::tibble(protein_id = "a1", sequence = "MKVLAAG"))
  ps2 <- new_protein_set("b", tibble::tibble(protein_id = "b1", sequence = "MKVLAAG"))
  pc <- new_protein_collection(list(ps1, ps2))
  expect_error(cluster_proteins(pc, timeout = -1), "timeout")
  expect_error(cluster_proteins(pc, timeout = "soon"), "timeout")
})

test_that("check_diamond_result reports timeouts and failures", {
  expect_error(
    check_diamond_result(list(timeout = TRUE, status = NA, stderr = ""), "makedb", 10),
    "did not finish within 10 seconds"
  )
  expect_error(
    check_diamond_result(list(timeout = FALSE, status = 1L, stderr = "boom"), "blastp", 10),
    "boom"
  )
  expect_true(check_diamond_result(list(timeout = FALSE, status = 0L, stderr = ""), "blastp", 10))
})

test_that("cluster_proteins aborts on a timed-out DIAMOND run", {
  diamond_path <- get_diamond_path()
  skip_if(is.null(diamond_path), "DIAMOND not installed")
  skip_on_os("windows")

  # A stand-in 'diamond' that never returns
  fake <- withr::local_tempfile(fileext = ".sh")
  writeLines(c("#!/bin/sh", "sleep 30"), fake)
  Sys.chmod(fake, "0755")

  p <- load_proteins(get_testdata_path(), quiet = TRUE)
  expect_error(
    cluster_proteins(p, tool_path = fake, timeout = 1, quiet = TRUE),
    "did not finish within"
  )
})

# Missing score files (issue 5) ------------------------------------------------

test_that("load_proteins warns when a score file is missing and still loads", {
  d <- copy_testdata()
  unlink(file.path(d, "assembly2_scored.csv"))

  expect_warning(
    p <- load_proteins(d, score_dir = d, quiet = TRUE),
    "assembly2"
  )
  expect_equal(unname(p$summary$has_scores), c(TRUE, FALSE, TRUE))
})

test_that("build_pa_matrix(score_threshold) works when one assembly has no scores", {
  ps1 <- new_protein_set("asm1", tibble::tibble(
    protein_id = c("asm1_p1", "asm1_p2"), sequence = c("MKV", "MKA"),
    custom_score = c(8, 2)
  ))
  ps2 <- new_protein_set("asm2", tibble::tibble(
    protein_id = c("asm2_p1", "asm2_p2"), sequence = c("MKV", "MKA")
  ))
  proteins <- new_protein_collection(list(ps1, ps2))
  og <- tibble::tibble(
    orthogroup_id = c("OG0001", "OG0001", "OG0002", "OG0002"),
    assembly = c("asm1", "asm2", "asm1", "asm2"),
    protein_id = c("asm1_p1", "asm2_p1", "asm1_p2", "asm2_p2")
  )
  singles <- tibble::tibble(assembly = "asm2", protein_id = "asm2_p3")
  ort <- new_orthogroup_result(og, method = "test", singletons = singles)

  expect_no_warning(pa <- build_pa_matrix(ort, proteins = proteins, score_threshold = 5))
  # OG0001 passes (8), OG0002 fails (2), the unscored singleton is dropped
  expect_equal(rownames(pa$matrix), "OG0001")
})

test_that("print.protein_collection reports partial scores", {
  ps1 <- new_protein_set("a", tibble::tibble(protein_id = "a1", sequence = "MKV", custom_score = 1))
  ps2 <- new_protein_set("b", tibble::tibble(protein_id = "b1", sequence = "MKV"))
  pc <- new_protein_collection(list(ps1, ps2))
  expect_output(print(pc), "Scores: present for 1 of 2 assemblies")
  expect_output(print(new_protein_collection(list(ps2))), "Scores: none loaded")
})

# Tool resolution (issues 6, 7) ------------------------------------------------

test_that("find_tool expands ~ in tool_path and conda_prefix", {
  home <- Sys.getenv("HOME")
  skip_if(!nzchar(home) || !dir.exists(home) || file.access(home, 2) != 0,
          "HOME not writable")

  tmp_dir <- file.path(home, paste0(".paneffectR-test-", Sys.getpid()))
  dir.create(file.path(tmp_dir, "bin"), recursive = TRUE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))
  tool <- file.path(tmp_dir, "bin", "fakediamond")
  writeLines("#!/bin/sh", tool)

  tilde_tool <- file.path("~", basename(tmp_dir), "bin", "fakediamond")
  tilde_prefix <- file.path("~", basename(tmp_dir))

  expect_false(grepl("~", find_tool("fakediamond", tool_path = tilde_tool), fixed = TRUE))
  expect_false(grepl("~", find_tool("fakediamond", conda_prefix = tilde_prefix), fixed = TRUE))
  expect_true(file.exists(find_tool("fakediamond", tool_path = tilde_tool)))
})

test_that("cluster_proteins accepts conda_env and reports when it is not found", {
  ps1 <- new_protein_set("a", tibble::tibble(protein_id = "a1", sequence = "MKVLAAG"))
  ps2 <- new_protein_set("b", tibble::tibble(protein_id = "b1", sequence = "MKVLAAG"))
  pc <- new_protein_collection(list(ps1, ps2))
  expect_error(
    cluster_proteins(pc, conda_env = "paneffectR-no-such-env-xyz"),
    "not found in conda env"
  )
})

# Score CSV validation (issue 8) -----------------------------------------------

test_that("load_scores aborts when required columns are missing", {
  tf <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("protein_id,other", "p1,1"), tf)
  expect_error(load_scores(tf), "custom_score")

  tf2 <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("id,custom_score", "p1,1"), tf2)
  expect_error(load_scores(tf2), "protein_id")
})

test_that("load_scores accepts a minimal file without readr parser warnings", {
  tf <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("protein_id,custom_score", "p1,1.5", "p2,RESULTS_MISSING"), tf)
  expect_no_warning(scores <- load_scores(tf))
  expect_equal(scores$custom_score, c(1.5, NA))
})

test_that("load_proteins warns when no score IDs match the FASTA", {
  d <- copy_testdata()
  csv <- file.path(d, "assembly1_scored.csv")
  scores <- readr::read_csv(csv, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  scores$protein_id <- paste0("other_", scores$protein_id)
  readr::write_csv(scores, csv)

  w <- capture_warnings(
    load_proteins(d, score_dir = d, pattern = "assembly1.faa", quiet = TRUE)
  )
  expect_true(any(grepl("No protein IDs", w)))
})

# Gzipped FASTA (issue 12) -----------------------------------------------------

test_that("gzipped FASTA files give the right assembly name and find score files", {
  d <- copy_testdata()
  for (f in list.files(d, pattern = "\\.faa$", full.names = TRUE)) {
    con <- gzfile(paste0(f, ".gz"), "w")
    writeLines(readLines(f), con)
    close(con)
    unlink(f)
  }

  p <- load_proteins(d, score_dir = d, pattern = "*.faa.gz", quiet = TRUE)
  expect_equal(p$summary$assembly_name, c("assembly1", "assembly2", "assembly3"))
  expect_true(all(p$summary$has_scores))
})

test_that("fasta_assembly_name strips compression and FASTA extensions", {
  expect_equal(fasta_assembly_name("/x/asm.faa"), "asm")
  expect_equal(fasta_assembly_name("/x/asm.faa.gz"), "asm")
  expect_equal(fasta_assembly_name("asm.v2.fasta.bz2"), "asm.v2")
})

# Empty / single-assembly input (issue 15) -------------------------------------

test_that("load_proteins warns on empty FASTA and single-assembly input", {
  d <- withr::local_tempdir()
  file.copy(get_testdata_path("assembly1.faa"), d)
  writeLines(character(), file.path(d, "empty.faa"))

  w <- character()
  withCallingHandlers(
    load_proteins(d, quiet = TRUE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("no sequences", w)))

  d2 <- withr::local_tempdir()
  file.copy(get_testdata_path("assembly1.faa"), d2)
  expect_warning(load_proteins(d2, quiet = TRUE), "Only one assembly")
})

test_that("load_proteins prints a summary unless quiet", {
  d <- copy_testdata()
  expect_message(load_proteins(d, score_dir = d), "Loaded 3 assemblies")
  expect_no_message(load_proteins(d, score_dir = d, quiet = TRUE))
})

# plot_heatmap performance (issue 10) ------------------------------------------

brute_jaccard <- function(mat) {
  n <- nrow(mat); d <- matrix(0, n, n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    a <- mat[i, ]; b <- mat[j, ]
    inter <- sum(a > 0 & b > 0); uni <- sum(a > 0 | b > 0)
    d[i, j] <- d[j, i] <- if (uni == 0) 0 else 1 - inter / uni
  }
  stats::as.dist(d)
}

brute_bray <- function(mat) {
  n <- nrow(mat); d <- matrix(0, n, n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    a <- mat[i, ]; b <- mat[j, ]
    den <- sum(a + b)
    d[i, j] <- d[j, i] <- if (den == 0) 0 else sum(abs(a - b)) / den
  }
  stats::as.dist(d)
}

test_that("vectorised Jaccard and Bray-Curtis match the brute-force definitions", {
  set.seed(42)
  m <- matrix(sample(0:3, 60, replace = TRUE), nrow = 15)
  m[3, ] <- 0  # all-zero row exercises the 0/0 guard
  expect_equal(as.vector(jaccard_dist(m)), as.vector(brute_jaccard(m)))
  expect_equal(as.vector(bray_curtis_dist(m)), as.vector(brute_bray(m)))
})

test_that("pattern_row_order groups identical rows together", {
  m <- rbind(c(1, 0), c(0, 1), c(1, 0), c(1, 1), c(0, 1))
  res <- pattern_row_order(m, "binary")
  expect_setequal(res$order, 1:5)
  # rows sharing a pattern are adjacent in the order
  ordered_patterns <- res$pattern[res$order]
  expect_equal(length(rle(ordered_patterns)$values), length(unique(res$pattern)))
})

test_that("plot_heatmap handles large matrices quickly via pattern clustering", {
  set.seed(1)
  n <- 20000
  m <- matrix(rbinom(n * 6, 1, 0.6), nrow = n,
              dimnames = list(sprintf("OG%05d", seq_len(n)), paste0("asm", 1:6)))
  storage.mode(m) <- "integer"
  pa <- new_pa_matrix(
    matrix = m,
    orthogroups = tibble::tibble(orthogroup_id = rownames(m), size = rowSums(m)),
    assemblies = tibble::tibble(assembly_name = colnames(m), n_orthogroups = colSums(m))
  )

  elapsed <- system.time(
    expect_message(ht <- plot_heatmap(pa), "distinct presence/absence pattern")
  )[["elapsed"]]
  expect_s4_class(ht, "Heatmap")
  expect_lt(elapsed, 30)

  # Below the threshold the full dendrogram path is used, silently
  expect_no_message(plot_heatmap(pa[1:50, ]))
})

# plot_upset warnings (issue 16) -----------------------------------------------

test_that("plot_upset does not emit UpSetR deprecation warnings", {
  visual_dir <- system.file("testdata", "visual", package = "paneffectR")
  if (visual_dir == "") visual_dir <- test_path("../../inst/testdata/visual")
  pa <- readRDS(file.path(visual_dir, "pa_binary.rds"))
  expect_no_warning(plot_upset(pa))
})
