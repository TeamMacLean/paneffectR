# Changelog

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
- Support for OrthoFinder and MMseqs2 backends (requires external
  installation)
- Import precomputed ortholog definitions

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
