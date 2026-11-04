# Load proteins from multiple assemblies

Main entry point for loading protein data. Discovers and loads FASTA
files, optionally with associated score and name mapping files.

## Usage

``` r
load_proteins(
  fasta_dir,
  score_dir = NULL,
  pattern = "*.faa",
  name_mapping = TRUE,
  quiet = FALSE
)
```

## Arguments

- fasta_dir:

  Character. Directory containing FASTA files.

- score_dir:

  Character. Optional directory containing score CSV files.

- pattern:

  Character. Glob pattern for FASTA files (default: "\*.faa"). Use e.g.
  `"*.faa.gz"` for gzipped files or `"*.fa*"` for mixed extensions.

- name_mapping:

  Logical. Whether to record name mapping files.

- quiet:

  Logical. If `FALSE` (default), print a one-line summary of what was
  loaded. Warnings are shown regardless.

## Value

A `protein_collection` object.

## File naming conventions

Each FASTA file in `fasta_dir` is treated as one assembly. The assembly
name is the file name without its extension (and without a trailing
`.gz`), so `isolate_A.faa` and `isolate_A.faa.gz` both give the assembly
`isolate_A`.

Optional companion files are matched by assembly name:

- `<assembly>_scored.csv` in `score_dir` - effector prediction scores
  (see below). A warning is issued for every assembly with no score
  file.

- `<assembly>_name_mapping.csv` in `fasta_dir` - omnieff name mapping
  (recorded in `source_files`, not merged).

## FASTA requirements

Files must contain protein (amino acid) sequences. Nucleotide input is
rejected with an error. Headers are split at the first whitespace: the
first token becomes `protein_id` and the remainder (if any) is kept in a
`description` column. For example `>XP_001 cytochrome P450 [P. oryzae]`
gives `protein_id = "XP_001"` and
`description = "cytochrome P450 [P. oryzae]"`. Protein IDs must be
unique across all assemblies. Trailing stop codons (`*`) are removed and
sequences are upper-cased.

## Score file format

A score file is a CSV with one row per protein and at least these
columns:

- `protein_id` - must match the FASTA IDs (first header token) exactly

- `custom_score` - numeric effector score used by `score_threshold` in
  [`build_pa_matrix()`](https://TeamMacLean.github.io/paneffectR/reference/build_pa_matrix.md)
  and by
  [`plot_scores()`](https://TeamMacLean.github.io/paneffectR/reference/plot_scores.md)

The omnieff pipeline also writes `assembly`, `score_rank`, and many
per-tool columns; all extra columns are kept. The value
`RESULTS_MISSING` (and empty strings) are read as `NA`. A warning is
issued if none of the IDs in a score file match the FASTA.

## Examples

``` r
# Bundled example data: three assemblies with scores
data_dir <- system.file("testdata", package = "paneffectR")
proteins <- load_proteins(fasta_dir = data_dir, score_dir = data_dir)
#> ✔ Loaded 3 assemblies (300 proteins) from
#>   /home/runner/work/_temp/Library/paneffectR/testdata
#> ℹ Scores: found for 3 of 3 assemblies, matched 300 of 300 proteins
proteins
#> -- protein_collection --
#> 3 assemblies, 300 total proteins
#> Scores: present
#> 
#> # A tibble: 3 × 3
#>   assembly_name n_proteins has_scores
#>   <chr>              <int> <lgl>     
#> 1 assembly1            100 TRUE      
#> 2 assembly2            100 TRUE      
#> 3 assembly3            100 TRUE      

if (FALSE) { # \dontrun{
# Load from omnieff output
proteins <- load_proteins(
  fasta_dir = "reformatted/",
  score_dir = "scored/"
)

# Load gzipped FASTAs without scores
proteins <- load_proteins(fasta_dir = "my_fastas/", pattern = "*.faa.gz")
} # }
```
