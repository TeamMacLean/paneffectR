# Load a single FASTA file

Parses a protein FASTA file into a tibble of protein IDs and sequences.
Headers are split at the first whitespace: the first token becomes
`protein_id` and the remainder is kept as `description` (`NA` if the
header had no description). Trailing stop codons (`*`) are removed and
sequences are upper-cased. Gzipped files are supported.

## Usage

``` r
load_fasta(path)
```

## Arguments

- path:

  Character. Path to FASTA file.

## Value

A tibble with columns: protein_id, sequence, description.

## Examples

``` r
fasta <- system.file("testdata", "assembly1.faa", package = "paneffectR")
seqs <- load_fasta(fasta)
head(seqs)
#> # A tibble: 6 × 3
#>   protein_id       sequence                                          description
#>   <chr>            <chr>                                             <chr>      
#> 1 assembly1_000039 MPFVLNLLVLVAVWYGLVRMAMLLKDDSEIKSFDPFAILGISAGVSER… NA         
#> 2 assembly1_000040 MEVAECTHKLASTTFVKTVEAPALAIMEPLLKRALAEGKTAVKRQASV… NA         
#> 3 assembly1_000094 MTDKDKPKDDDDCYVKVNTPGGVGDTNSPNKGGLPFTPSNLCIKNLEY… NA         
#> 4 assembly1_000122 MGSSPIFRPFKCKFSMFSKSILIFLLILPLISVILLSFSKNEKFIKNF… NA         
#> 5 assembly1_000164 MLFCVMVGTGCQLLGMALVTLFFAAVGVLAPSNRGKLVIALLVCFVLL… NA         
#> 6 assembly1_000194 MMGNLAKDIVNGHREKTLALLWKLISCFQLQALVDAQTMRREISNVVK… NA         
```
