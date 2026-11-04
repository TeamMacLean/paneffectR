# Changelog

## paneffectR 0.1.2

### Heatmap styles

- [`plot_heatmap()`](https://TeamMacLean.github.io/paneffectR/reference/plot_heatmap.md)
  gains a `style` argument with three ready-made looks: `"clean"` (new
  default: one hue on a light grey grid), `"pan"` (rows split into Core
  / Accessory / Unique blocks with a coloured side bar; binary cells
  take the class colour) and `"mono"` (black on white for print).
- All styles label the legend “Absent”/“Present” with a meaningful
  title, draw a light grid between cells on matrices of up to 200 rows,
  use slimmer dendrograms and angled column labels. Score matrices use
  an orange ramp so they are not mistaken for presence/absence; count
  matrices whose maximum is 1 use the binary colours.
- Arguments passed through `...` to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  now override the style’s settings instead of erroring on duplicate
  arguments.
- New vignette, “Styling Heatmaps”, covering the three styles, shared
  arguments, saving figures and how to go further with ComplexHeatmap.
  README and the other vignettes use the new styles.

## paneffectR 0.1.1

Fixes from a new-user functional test pass (see GitHub issues for
details).

### Bug fixes

- FASTA headers with descriptions (e.g. NCBI, UniProt, BRAKER output) no
  longer silently break clustering.
  [`load_fasta()`](https://TeamMacLean.github.io/paneffectR/reference/load_fasta.md)
  now uses the first whitespace-delimited token as `protein_id` and
  keeps the remainder in a new `description` column.
  [`cluster_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/cluster_proteins.md)
  also validates that DIAMOND’s IDs map back to loaded proteins.
- Nucleotide FASTA input (which made DIAMOND hang indefinitely) is
  detected and rejected with a clear error, both at load time and before
  clustering. External tool calls now have a `timeout` (default 3600 s).
- `build_pa_matrix(score_threshold = )` no longer crashes when some
  assemblies have no score file; orthogroups without any score are
  dropped rather than erroring.
- `tool_path` and `conda_prefix` containing `~` are expanded before use.
- [`cluster_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/cluster_proteins.md)
  gains `conda_env` to use a named conda/mamba environment, matching the
  documented resolution order.
- Gzipped FASTA files (`*.faa.gz`) now yield the correct assembly name
  so score files are found.
- Trailing stop codons (`*`) are stripped and sequences upper-cased on
  load.
- [`plot_upset()`](https://TeamMacLean.github.io/paneffectR/reference/plot_upset.md)
  no longer emits ggplot2 deprecation warnings from UpSetR.

### Improvements

- [`load_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/load_proteins.md)
  prints a one-line summary and warns about missing score files, score
  IDs that match no FASTA IDs, empty FASTA files and single-assembly
  input.
- [`load_scores()`](https://TeamMacLean.github.io/paneffectR/reference/load_scores.md)
  errors clearly when `protein_id` or `custom_score` is missing.
- [`cluster_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/cluster_proteins.md)
  reports progress for each step and the temporary directory when
  `keep_temp = TRUE`. Use `quiet = TRUE` to silence.
- [`plot_heatmap()`](https://TeamMacLean.github.io/paneffectR/reference/plot_heatmap.md)
  is much faster: Jaccard and Bray-Curtis distances are vectorised, and
  matrices with more than `max_cluster_rows` rows (default 5000) are
  ordered by clustering distinct presence/absence patterns.
- [`print()`](https://rdrr.io/r/base/print.html) for
  `protein_collection` reports whether scores are loaded.

### Documentation

- README install instructions use `pak`/`remotes` (devtools \>= 2.5 no
  longer bundles `install_github()`), and
  [`build_pa_matrix()`](https://TeamMacLean.github.io/paneffectR/reference/build_pa_matrix.md)
  examples pass `proteins`.
- [`?load_proteins`](https://TeamMacLean.github.io/paneffectR/reference/load_proteins.md)
  documents file naming conventions and the score CSV format.
- All exported functions now have runnable examples.
- Fixed non-running vignette code (`row.names` typo; hard-coded assembly
  name).

### Breaking changes

- `filter_by_score()` (an unimplemented stub) has been removed. Use
  `build_pa_matrix(score_threshold = )`.
- `cluster_proteins(method = )` only accepts `"diamond_rbh"`; the
  unimplemented `"orthofinder"` and `"mmseqs2"` values are no longer
  listed.
- [`load_fasta()`](https://TeamMacLean.github.io/paneffectR/reference/load_fasta.md)
  returns an additional `description` column.

## paneffectR 0.1.0

Initial release of paneffectR. \## Features

#### Data Loading

- [`load_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/load_proteins.md)
  loads multiple FASTA files with optional effector scores
- [`load_fasta()`](https://TeamMacLean.github.io/paneffectR/reference/load_fasta.md),
  [`load_scores()`](https://TeamMacLean.github.io/paneffectR/reference/load_scores.md),
  [`load_name_mapping()`](https://TeamMacLean.github.io/paneffectR/reference/load_name_mapping.md)
  for individual files
- Support for omnieff pipeline output format

#### Clustering

- [`cluster_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/cluster_proteins.md)
  with DIAMOND RBH + expansion algorithm
- Configurable identity, coverage, and e-value thresholds

#### Matrix Building

- [`build_pa_matrix()`](https://TeamMacLean.github.io/paneffectR/reference/build_pa_matrix.md)
  creates presence/absence matrices from clustering results
- Binary, count, and score-based matrix types
- Score threshold filtering for effector analysis
- Optional singleton inclusion

#### Visualization

- [`plot_heatmap()`](https://TeamMacLean.github.io/paneffectR/reference/plot_heatmap.md)
  with ComplexHeatmap backend
  - Row and column clustering
  - Multiple distance metrics (binary, Jaccard, Bray-Curtis)
  - Row and column annotations
- [`plot_upset()`](https://TeamMacLean.github.io/paneffectR/reference/plot_upset.md)
  for set intersection visualization
- [`plot_scores()`](https://TeamMacLean.github.io/paneffectR/reference/plot_scores.md)
  for effector score distributions
- [`plot_dendro()`](https://TeamMacLean.github.io/paneffectR/reference/plot_dendro.md)
  for assembly clustering dendrograms

#### Utilities

- [`get_singletons()`](https://TeamMacLean.github.io/paneffectR/reference/get_singletons.md),
  [`n_singletons()`](https://TeamMacLean.github.io/paneffectR/reference/n_singletons.md),
  [`singletons_by_assembly()`](https://TeamMacLean.github.io/paneffectR/reference/singletons_by_assembly.md)
  for singleton analysis
- [`check_tool_installed()`](https://TeamMacLean.github.io/paneffectR/reference/check_tool_installed.md),
  [`find_tool()`](https://TeamMacLean.github.io/paneffectR/reference/find_tool.md)
  for external tool management

### Documentation

- Getting Started vignette with core workflow
- Effector Analysis vignette for omnieff users
- Pan-Genome Analysis vignette for general comparisons
- Algorithm Deep Dive vignette with technical details
