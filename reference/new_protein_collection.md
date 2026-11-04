# Create a protein_collection object

Container for multiple assemblies.

## Usage

``` r
new_protein_collection(assemblies)
```

## Arguments

- assemblies:

  List of protein_set objects.

## Value

A `protein_collection` S3 object.

## Examples

``` r
ps_a <- new_protein_set("isolate_A",
  tibble::tibble(protein_id = "A_001", sequence = "MKVLAAGIV"))
ps_b <- new_protein_set("isolate_B",
  tibble::tibble(protein_id = "B_001", sequence = "MKVLAAGLV"))
pc <- new_protein_collection(list(ps_a, ps_b))
pc
#> -- protein_collection --
#> 2 assemblies, 2 total proteins
#> Scores: none loaded
#> 
#> # A tibble: 2 × 3
#>   assembly_name n_proteins has_scores
#>   <chr>              <int> <lgl>     
#> 1 isolate_A              1 FALSE     
#> 2 isolate_B              1 FALSE     
```
