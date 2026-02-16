# Package index

## Loading Data

Functions for loading protein sequences and effector prediction scores
from FASTA files and omnieff output.

- [`load_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/load_proteins.md)
  : Load proteins from multiple assemblies
- [`load_fasta()`](https://TeamMacLean.github.io/paneffectR/reference/load_fasta.md)
  : Load a single FASTA file
- [`load_scores()`](https://TeamMacLean.github.io/paneffectR/reference/load_scores.md)
  : Load effector scores from CSV
- [`load_name_mapping()`](https://TeamMacLean.github.io/paneffectR/reference/load_name_mapping.md)
  : Load name mapping from CSV

## Creating Objects

Constructors for the core S3 classes used throughout the package.

- [`new_protein_set()`](https://TeamMacLean.github.io/paneffectR/reference/new_protein_set.md)
  : Create a protein_set object
- [`new_protein_collection()`](https://TeamMacLean.github.io/paneffectR/reference/new_protein_collection.md)
  : Create a protein_collection object
- [`new_orthogroup_result()`](https://TeamMacLean.github.io/paneffectR/reference/new_orthogroup_result.md)
  : Create an orthogroup_result object
- [`new_pa_matrix()`](https://TeamMacLean.github.io/paneffectR/reference/new_pa_matrix.md)
  : Create a pa_matrix object

## Clustering

Functions for clustering proteins into orthogroups using sequence
similarity.

- [`cluster_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/cluster_proteins.md)
  : Cluster proteins across assemblies

## Matrix Building

Build presence/absence matrices from clustering results.

- [`build_pa_matrix()`](https://TeamMacLean.github.io/paneffectR/reference/build_pa_matrix.md)
  : Build presence/absence matrix from clustering results
- [`filter_by_score()`](https://TeamMacLean.github.io/paneffectR/reference/filter_by_score.md)
  : Filter presence/absence matrix by score
- [`` `[`( ``*`<pa_matrix>`*`)`](https://TeamMacLean.github.io/paneffectR/reference/sub-.pa_matrix.md)
  : Subset a pa_matrix object
- [`as.data.frame(`*`<pa_matrix>`*`)`](https://TeamMacLean.github.io/paneffectR/reference/as.data.frame.pa_matrix.md)
  : Convert pa_matrix to data.frame

## Singletons

Functions for working with unclustered proteins (singletons).

- [`get_singletons()`](https://TeamMacLean.github.io/paneffectR/reference/get_singletons.md)
  : Get singletons from clustering result
- [`n_singletons()`](https://TeamMacLean.github.io/paneffectR/reference/n_singletons.md)
  : Count singletons
- [`singletons_by_assembly()`](https://TeamMacLean.github.io/paneffectR/reference/singletons_by_assembly.md)
  : Summarize singletons by assembly

## Visualization

Publication-ready plots for exploring presence/absence patterns.

- [`plot_heatmap()`](https://TeamMacLean.github.io/paneffectR/reference/plot_heatmap.md)
  : Plot presence/absence heatmap
- [`plot_upset()`](https://TeamMacLean.github.io/paneffectR/reference/plot_upset.md)
  : Plot UpSet diagram of orthogroup sharing
- [`plot_scores()`](https://TeamMacLean.github.io/paneffectR/reference/plot_scores.md)
  : Plot effector score distributions
- [`plot_dendro()`](https://TeamMacLean.github.io/paneffectR/reference/plot_dendro.md)
  : Plot assembly clustering dendrogram
- [`plot_pan_structure()`](https://TeamMacLean.github.io/paneffectR/reference/plot_pan_structure.md)
  : Plot pan-genome structure
- [`plot_assembly_composition()`](https://TeamMacLean.github.io/paneffectR/reference/plot_assembly_composition.md)
  : Plot assembly composition

## Utilities

Helper functions for external tool management.

- [`check_tool_installed()`](https://TeamMacLean.github.io/paneffectR/reference/check_tool_installed.md)
  : Check if an external tool is installed
- [`find_tool()`](https://TeamMacLean.github.io/paneffectR/reference/find_tool.md)
  : Find an external tool
