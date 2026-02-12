# Clustering Verification Test Data

This directory contains real protein sequences from well-annotated fungal species with known orthology relationships. It is used to verify that the paneffectR clustering pipeline produces biologically correct results.

## Overview

- **3 assemblies**: S. cerevisiae (scer), S. pombe (spom), N. crassa (ncra)
- **16 proteins total**: 5 conserved protein families + 1 species-specific singleton
- **5 expected orthogroups**: Each conserved protein family should cluster together
- **1 singleton**: Species-specific protein should not cluster

## Organisms

| Code | Species | Common Name | Why Selected |
|------|---------|-------------|--------------|
| scer | *Saccharomyces cerevisiae* (S288c) | Budding yeast | Gold standard genome annotation |
| spom | *Schizosaccharomyces pombe* (972) | Fission yeast | Well-annotated, diverged ~400 MYA from S. cerevisiae |
| ncra | *Neurospora crassa* (74-OR23-1A) | Red bread mold | Filamentous fungus, diverged ~800 MYA from yeasts |

## Protein Families

### Conserved Orthologs (should cluster together)

| Expected OG | Function | Conservation |
|-------------|----------|--------------|
| OG_actin | Cytoskeleton | Extremely high (>85% identity across fungi) |
| OG_gapdh | Glycolysis | Very high (~70% identity) |
| OG_histone_h3 | Chromatin | Extremely high (>90% identity) |
| OG_ef1a | Translation | Very high (~80% identity) |
| OG_cyc | Electron transport | High (~60% identity) |

### Singleton (should not cluster)

| Protein | Function | Notes |
|---------|----------|-------|
| scer_MFA1 | Mating pheromone | Species-specific, no orthologs in fission yeast or N. crassa |

## Sequence Sources

All sequences were retrieved from UniProt (Swiss-Prot reviewed entries) on 2026-02-12.

| Protein ID | UniProt | Organism | Gene |
|------------|---------|----------|------|
| scer_ACT1 | [P60010](https://www.uniprot.org/uniprotkb/P60010) | S. cerevisiae | ACT1 |
| spom_act1 | [P10989](https://www.uniprot.org/uniprotkb/P10989) | S. pombe | act1 |
| ncra_act | [P78711](https://www.uniprot.org/uniprotkb/P78711) | N. crassa | act |
| scer_TDH1 | [P00360](https://www.uniprot.org/uniprotkb/P00360) | S. cerevisiae | TDH1 |
| spom_tdh1 | [P78958](https://www.uniprot.org/uniprotkb/P78958) | S. pombe | tdh1 |
| ncra_gpd1 | [P54118](https://www.uniprot.org/uniprotkb/P54118) | N. crassa | gpd-1 |
| scer_HHT1 | [P61830](https://www.uniprot.org/uniprotkb/P61830) | S. cerevisiae | HHT1 |
| spom_hht1 | [P09988](https://www.uniprot.org/uniprotkb/P09988) | S. pombe | hht1 |
| ncra_hh3 | [P07041](https://www.uniprot.org/uniprotkb/P07041) | N. crassa | hh3 |
| scer_TEF1 | [P02994](https://www.uniprot.org/uniprotkb/P02994) | S. cerevisiae | TEF1 |
| spom_tef101 | [P0CT53](https://www.uniprot.org/uniprotkb/P0CT53) | S. pombe | tef101 |
| ncra_tef1 | [Q01372](https://www.uniprot.org/uniprotkb/Q01372) | N. crassa | tef-1 |
| scer_CYC1 | [P00044](https://www.uniprot.org/uniprotkb/P00044) | S. cerevisiae | CYC1 |
| spom_cyc1 | [P00046](https://www.uniprot.org/uniprotkb/P00046) | S. pombe | cyc1 |
| ncra_cyc1 | [P00048](https://www.uniprot.org/uniprotkb/P00048) | N. crassa | cyc-1 |
| scer_MFA1 | [P01149](https://www.uniprot.org/uniprotkb/P01149) | S. cerevisiae | MF(ALPHA)1 |

## Files

| File | Description |
|------|-------------|
| `scer.faa` | S. cerevisiae protein sequences (6 proteins) |
| `spom.faa` | S. pombe protein sequences (5 proteins) |
| `ncra.faa` | N. crassa protein sequences (5 proteins) |
| `expected_orthogroups.csv` | Expected clustering results with protein metadata |
| `known_orthologs.rds` | Pre-built protein_collection for testing |
| `README.md` | This file |

## Expected Results

When running the clustering pipeline:

```r
proteins <- readRDS("inst/testdata/clustering/known_orthologs.rds")
clusters <- cluster_proteins(proteins, method = "diamond_rbh")
```

The expected results are:

1. **5 orthogroups** should be detected (actin, GAPDH, histone H3, EF-1α, cytochrome c)
2. **Each orthogroup** should contain exactly 3 proteins (one per species)
3. **scer_MFA1** should be a singleton (not assigned to any orthogroup)

## Verification Checklist

- [ ] All clear orthologs cluster into the same orthogroup
- [ ] Each orthogroup has exactly 3 members (one per species)
- [ ] scer_MFA1 is identified as a singleton
- [ ] PA matrix shows expected presence pattern:
  - Conserved proteins: present in all 3 assemblies
  - Singleton: present only in scer

## Caveats

1. **Paralogs**: S. cerevisiae has multiple GAPDH genes (TDH1, TDH2, TDH3) and histone H3 genes (HHT1, HHT2). We only include one representative from each family.

2. **Gene names**: Gene names differ between species (e.g., TDH1 in yeast, gpd-1 in N. crassa). This is expected for orthologs.

3. **Sequence divergence**: While these proteins are conserved, there is still ~10-40% sequence divergence. The clustering algorithm must be sensitive enough to detect homology across this divergence.

## Regenerating Test Data

If sequences need to be updated:

1. Search UniProt for each protein using organism and gene name
2. Download FASTA sequences from Swiss-Prot (reviewed entries)
3. Rename headers to use assembly prefix (e.g., `scer_ACT1`)
4. Update `expected_orthogroups.csv` with new accessions if changed
5. Rebuild `known_orthologs.rds` using:

```r
devtools::load_all()
proteins <- load_proteins(
  fasta_dir = "inst/testdata/clustering",
  pattern = "*.faa"
)
saveRDS(proteins, "inst/testdata/clustering/known_orthologs.rds")
```
