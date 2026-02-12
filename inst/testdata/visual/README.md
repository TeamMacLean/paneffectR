# Visual Test Data

Synthetic data designed for human verification of plots.

## Design

- **6 assemblies**: asm_A, asm_B, asm_C, asm_D, asm_E, asm_F
- **50 orthogroups** with predictable patterns:

| Category | Count | Present in | Score range | Expected appearance |
|----------|-------|------------|-------------|---------------------|
| core | 10 | All 6 | 7-10 | Solid row, hot color |
| accessory | 15 | 3-5 | 4-7 | Partial row, medium |
| rare | 15 | 2 | 2-4 | Sparse row, cool |
| unique | 10 | 1 only | 1-3 | Single cell per row |

## Verification Checklist

When reviewing plots:

### Heatmap
- [ ] Core orthogroups (OG_core_*) show solid horizontal bands
- [ ] Unique orthogroups show single cells
- [ ] Score heatmap shows color gradient (core=hot, unique=cool)

### UpSet
- [ ] Largest intersection is "all 6 assemblies" (10 orthogroups)
- [ ] Each assembly has 1-2 unique orthogroups
- [ ] Total adds up to 50

### Dendrogram
- [ ] Assemblies cluster by shared content
- [ ] asm_A and asm_B should cluster (share most accessory)

### Score Distribution
- [ ] Higher scores concentrated in some assemblies
- [ ] Threshold line visible when specified

## Files

- `proteins_visual.rds` - protein_collection object
- `clusters_visual.rds` - orthogroup_result object
- `pa_visual.rds` - pa_matrix object (pre-built for convenience)
