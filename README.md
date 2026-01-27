# paneffectR

Comparative genomics of effector proteins across multiple genome assemblies.

## Overview

**paneffectR** clusters proteins from different assemblies into ortholog groups, builds presence/absence matrices, and generates publication-ready visualizations. It is designed for effector protein analysis from the [omnieff](https://github.com/username/omnieff) pipeline but works with any protein sets.

## Installation

paneffectR requires Bioconductor packages. Install BiocManager first:
```r
install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "Biostrings"))
```

Then install paneffectR from GitHub:

```r
# install.packages("devtools")
devtools::install_github("username/paneffectR")
```

## External Dependencies

At least one clustering tool must be installed:

| Tool | Purpose | Installation |
|------|---------|--------------|
| DIAMOND | Fast protein alignment (default) | `mamba install -c bioconda diamond` |
| OrthoFinder | Comprehensive ortholog inference | `mamba install -c bioconda orthofinder` |
| MMseqs2 | Ultra-fast clustering | `mamba install -c bioconda mmseqs2` |

## Quick Start

```r
library(paneffectR)

# Load proteins from multiple assemblies
proteins <- load_proteins(
  fasta_dir = "path/to/fastas/",
  score_dir = "path/to/scores/"  # optional
)

# Cluster into orthogroups
clusters <- cluster_proteins(proteins, method = "diamond_rbh")

# Build presence/absence matrix
pa <- build_pa_matrix(clusters, score_threshold = 5.0)

# Visualize
plot_heatmap(pa)
plot_upset(pa)
```

## Use Cases

1. **Effector analysis**: Take omnieff output, find equivalent proteins across assemblies, filter by effector score
2. **Pan-genome analysis**: Compare any proteins across assemblies, identify core/accessory/unique proteins

## License

MIT
