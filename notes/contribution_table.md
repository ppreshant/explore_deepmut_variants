# Mapping contribution to each nucleotide position

## Goal
To map the additive contribution of each nucleotide position to the overall catalysis and leak scores (and net catalysis).

## Style
- Follow the tidyverse style and the existing code written in exloratory_analysis.qmd. 
  - Use dplyr and tibble functions where possible and vectorize where elegant for performance. 
  - Avoid for loops where possible.
- Write code in reusable functions where relevant/useful. 


## Plan 
1. Make a tibble with each row representing a nucleotide position and the corresponding nucleotide (A, C, G, T) at that position in the WT sequence (nham == 0). 
  - Create 4 columns for each nucleotide (A, C, G, T) and initialize them to 0. These columns will be used to count the number of variants where that WT nucleotide is mutated to a different nucleotide.
  - create 1 more column to represent the total number of variants where that WT nucleotide is mutated to a different nucleotide. Initialize it to 0.
  - Create additional columns for the sum of relative_catalysis, relative_leak, and net_catalysis for each nucleotide at that position. Initialize them to NA.
2. For each variant, identify the positions of the mutations (snps) and the corresponding nucleotides at those positions; Populate the counts and score columns for each nucleotide at each position in the tibble.
  - For the scores, subtract 1 before adding to the sum score columns, since the scores are normalized to the WT (which has a score of 1). This will allow us to calculate the average contribution of each nucleotide at each position to the overall catalysis and leak scores. 

## Useful pseudocode
This uses data.frame, but remake it using tibble and dplyr functions. 

```r
library(data.table)

wt <- toupper("YOUR_WT_SEQUENCE_HERE")   # same string you passed to --wildtypeSequence, uppercase,
                                          # and with constant regions removed if you used lowercase masking
wt_chars <- strsplit(wt, "")[[1]]

get_muts <- function(seq) {
  s <- strsplit(seq, "")[[1]]
  idx <- which(s != wt_chars)
  if (length(idx) == 0) return(NA_character_)
  paste0(wt_chars[idx], idx, s[idx], collapse = ";")   # e.g. "A12G;C15T"
}

dt <- as.data.table(all_variants)     # whatever your object is called
dt[, mutations := vapply(nt_seq, get_muts, character(1))]
dt[, positions := vapply(nt_seq, function(x) {
  s <- strsplit(x, "")[[1]]; paste(which(s != wt_chars), collapse = ";")
}, character(1))]
```