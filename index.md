# paneffectR

Comparative genomics of effector proteins across multiple genome
assemblies.

## Why paneffectR?

When studying plant pathogens, a key question is: **which effector
proteins are shared across isolates, and which are unique?**
Understanding effector repertoires helps identify core virulence
factors, track pathogen evolution, and discover candidate avirulence
genes.

**paneffectR** solves this by:

1.  **Clustering proteins into orthogroups** - Finding equivalent
    proteins across assemblies using sequence similarity
2.  **Building presence/absence matrices** - Creating structured data
    showing which proteins exist in which assemblies
3.  **Filtering by effector scores** - Focusing on high-confidence
    effector predictions (when using
    [omnieff](https://github.com/TeamMacLean/omnieff) output)
4.  **Generating publication-ready visualizations** - Heatmaps, UpSet
    plots, and dendrograms

While designed for effector analysis, paneffectR works with any protein
sets for general pan-genome comparisons.

## Installation

### R Dependencies

paneffectR requires Bioconductor packages. Install them first:

``` r
install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "Biostrings"))
```

Then install paneffectR from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("TeamMacLean/paneffectR")
```

### External Dependencies

At least one clustering tool must be installed. DIAMOND is recommended:

| Tool                                                    | Purpose                                  | Installation                            |
|---------------------------------------------------------|------------------------------------------|-----------------------------------------|
| [DIAMOND](https://github.com/bbuchfink/diamond)         | Fast protein alignment (default)         | `mamba install -c bioconda diamond`     |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | Comprehensive ortholog inference         | `mamba install -c bioconda orthofinder` |
| [MMseqs2](https://github.com/soedinglab/MMseqs2)        | Ultra-fast clustering for large datasets | `mamba install -c bioconda mmseqs2`     |

## Quick Start

``` r
library(paneffectR)

# Load proteins from multiple assemblies
proteins <- load_proteins(

  fasta_dir = "path/to/fastas/",
  score_dir = "path/to/scores/"
)

# Cluster into orthogroups
clusters <- cluster_proteins(proteins, method = "diamond_rbh")

# Build presence/absence matrix (filter to high-scoring effectors)
pa <- build_pa_matrix(clusters, score_threshold = 5.0)

# Visualize
ht <- plot_heatmap(pa)
ComplexHeatmap::draw(ht)

plot_upset(pa, min_size = 2)
```

## Use Cases

### Effector Comparative Genomics

Take output from the [omnieff](https://github.com/TeamMacLean/omnieff)
pipeline, find orthologous effectors across assemblies, and filter by
prediction confidence:

``` r
# Load omnieff output (FASTAs + scores)
proteins <- load_proteins(
  fasta_dir = "omnieff_output/reformatted/",
  score_dir = "omnieff_output/scored/"
)

# Cluster and build matrix with score threshold
clusters <- cluster_proteins(proteins)
pa <- build_pa_matrix(clusters, score_threshold = 5.0)

# Visualize effector repertoires
plot_heatmap(pa) |> ComplexHeatmap::draw()
```

### General Pan-Genome Analysis

Compare any protein sets without effector scores:

``` r
# Load raw FASTAs
proteins <- load_proteins(fasta_dir = "my_assemblies/")

# Binary presence/absence analysis
clusters <- cluster_proteins(proteins)
pa <- build_pa_matrix(clusters, type = "binary")

# Identify core vs accessory proteins
plot_upset(pa, min_size = 2)
plot_dendro(pa, distance_method = "jaccard")
```

## Documentation

- **[Getting
  Started](https://TeamMacLean.github.io/paneffectR/articles/getting-started.html)** -
  Core workflow tutorial
- **[Effector
  Analysis](https://TeamMacLean.github.io/paneffectR/articles/effector-analysis.html)** -
  Working with omnieff output
- **[Pan-Genome
  Analysis](https://TeamMacLean.github.io/paneffectR/articles/pan-genome.html)** -
  General protein comparisons
- **[Algorithm Deep
  Dive](https://TeamMacLean.github.io/paneffectR/articles/algorithms.html)** -
  Technical details for bioinformaticians
- **[Function
  Reference](https://TeamMacLean.github.io/paneffectR/reference/index.html)** -
  Complete API documentation

## Citation

If you use paneffectR in your research, please cite:

> MacLean, D. (2026). paneffectR: Comparative Genomics of Effector
> Proteins. R package version 0.1.0.
> <https://github.com/TeamMacLean/paneffectR>

## License

MIT
