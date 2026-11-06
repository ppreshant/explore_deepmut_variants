## Contribution of substitutions by nucleotide position ----

#' Summarise substitution counts and WT-normalised score changes by position.
#'
#' Each substitution in a multi-mutant receives the complete variant score
#' change. Indels are excluded because their positions cannot be aligned to WT.
make_contribution_table <- function(df) {

  # Identify the single WT sequence used as the score normalisation reference.
  wt_rows <- df |>
    filter(!is.na(Nham_nt), Nham_nt == 0) |>
    distinct(nt_seq)

  if (nrow(wt_rows) != 1) {
    stop("Expected exactly one unique wild-type sequence (Nham_nt == 0).")
  }

  wt_sequence <- toupper(wt_rows$nt_seq[[1]])
  wt_nucleotides <- strsplit(wt_sequence, "", fixed = TRUE)[[1]]

  # Keep aligned variants with complete scores and expand each mutation
  # into one row so multi-mutants contribute to every affected position.
  variants <-
    filter(df, !is.na(Nham_nt)) |> # remove indels (retain SNPs)
    mutate(
      .positions = map(nt_seq, \(sequence) {
        sequence_nucleotides <- strsplit(toupper(sequence), "", fixed = TRUE)[[1]]
        which(sequence_nucleotides != wt_nucleotides) # Locate substitutions.
      })
    ) |>
    unnest_longer(.positions, values_to = "position") |>
    mutate(
      wt_nucleotide = wt_nucleotides[position], # Record the WT base.
      mutant_nucleotide = str_sub(toupper(nt_seq), position, position),
      catalysis_change = relative_catalysis - 1, # Convert to WT-relative change.
      leak_change = relative_leak - 1, # Convert to WT-relative change.
      net_catalysis_change = net_catalysis - 1 # Convert to WT-relative change.
    ) |>
    filter(mutant_nucleotide %in% c("A", "C", "G", "T")) # Keep SNPs.

  # Aggregate mutation counts and score changes for each position/base pair.
  by_nucleotide <- variants |>
    group_by(position, wt_nucleotide, mutant_nucleotide) |>
    summarise(
      count = n(), # Number of variants with this substitution.
      relative_catalysis_sum = sum(catalysis_change), # Total catalysis change.
      relative_leak_sum = sum(leak_change), # Total leak change.
      net_catalysis_sum = sum(net_catalysis_change), # Total net change.
      relative_catalysis_mean = mean(catalysis_change), # Mean catalysis change.
      relative_leak_mean = mean(leak_change), # Mean leak change.
      net_catalysis_mean = mean(net_catalysis_change), # Mean net change.
      .groups = "drop"
    ) |>
    # Add zero-count rows so every position has all four destination bases.
    complete(
      position = seq_along(wt_nucleotides),
      mutant_nucleotide = c("A", "C", "G", "T"),
      fill = list(
        count = 0,
        relative_catalysis_sum = 0,
        relative_leak_sum = 0,
        net_catalysis_sum = 0,
        relative_catalysis_mean = NA_real_,
        relative_leak_mean = NA_real_,
        net_catalysis_mean = NA_real_
      )
    ) |>
    mutate(wt_nucleotide = wt_nucleotides[position]) |>
    # Put each destination base's metrics into separate columns.
    pivot_wider(
      names_from = mutant_nucleotide,
      values_from = c(
        count,
        relative_catalysis_sum,
        relative_leak_sum,
        net_catalysis_sum,
        relative_catalysis_mean,
        relative_leak_mean,
        net_catalysis_mean
      ),
      names_glue = "x>{mutant_nucleotide}_{.value}"
    )

  by_nucleotide |>
    mutate(
      total_mutations = rowSums(
        pick(matches("^x>[ACGT]_count$")),
        na.rm = TRUE
      ), # Total substitutions at this position.
      total_relative_catalysis_sum = rowSums(
        pick(matches("^x>[ACGT]_relative_catalysis_sum$")),
        na.rm = TRUE
      ), # Total catalysis change at this position.
      total_relative_leak_sum = rowSums(
        pick(matches("^x>[ACGT]_relative_leak_sum$")),
        na.rm = TRUE
      ), # Total leak change at this position.
      total_net_catalysis_sum = rowSums(
        pick(matches("^x>[ACGT]_net_catalysis_sum$")),
        na.rm = TRUE
      ), # Total net catalysis change at this position.
      total_relative_catalysis_mean = if_else(
        total_mutations > 0,
        total_relative_catalysis_sum / total_mutations,
        NA_real_
      ), # Mean catalysis change across all substitutions.
      total_relative_leak_mean = if_else(
        total_mutations > 0,
        total_relative_leak_sum / total_mutations,
        NA_real_
      ), # Mean leak change across all substitutions.
      total_net_catalysis_mean = if_else(
        total_mutations > 0,
        total_net_catalysis_sum / total_mutations,
        NA_real_
      ) # Mean net catalysis change across all substitutions.
    ) |>
    ungroup() |>
    select(
      position,
      wt_nucleotide,
      total_mutations,
      total_relative_catalysis_mean,
      total_relative_leak_mean,
      total_net_catalysis_mean,
      total_relative_catalysis_sum,
      total_relative_leak_sum,
      total_net_catalysis_sum,
      everything()
    )
}
