# paneffectR 0.1.0

Initial release of paneffectR.
## Features

### Data Loading
* `load_proteins()` loads multiple FASTA files with optional effector scores
* `load_fasta()`, `load_scores()`, `load_name_mapping()` for individual files
* Support for omnieff pipeline output format

### Clustering
* `cluster_proteins()` with DIAMOND RBH + expansion algorithm
* Configurable identity, coverage, and e-value thresholds
* Support for OrthoFinder and MMseqs2 backends (requires external installation)
* Import precomputed ortholog definitions

### Matrix Building
* `build_pa_matrix()` creates presence/absence matrices from clustering results
* Binary, count, and score-based matrix types
* Score threshold filtering for effector analysis
* Optional singleton inclusion

### Visualization
* `plot_heatmap()` with ComplexHeatmap backend
  - Row and column clustering
  - Multiple distance metrics (binary, Jaccard, Bray-Curtis)
  - Row and column annotations
* `plot_upset()` for set intersection visualization
* `plot_scores()` for effector score distributions
* `plot_dendro()` for assembly clustering dendrograms

### Utilities
* `get_singletons()`, `n_singletons()`, `singletons_by_assembly()` for singleton analysis
* `check_tool_installed()`, `find_tool()` for external tool management

## Documentation

* Getting Started vignette with core workflow
* Effector Analysis vignette for omnieff users
* Pan-Genome Analysis vignette for general comparisons
* Algorithm Deep Dive vignette with technical details
