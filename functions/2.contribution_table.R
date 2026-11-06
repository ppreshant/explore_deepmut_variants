## Contribution of substitutions by nucleotide position ----

#' Summarise substitution counts and WT-normalised score changes by position.
#'
#' Each substitution in a multi-mutant receives the complete variant score
#' change. Indels are excluded because their positions cannot be aligned to WT.
make_contribution_table <- function(df) {

  wt_rows <- df |>
    dplyr::filter(!is.na(Nham_nt), Nham_nt == 0) |>
    dplyr::distinct(nt_seq)

  if (nrow(wt_rows) != 1) {
    stop("Expected exactly one unique wild-type sequence (Nham_nt == 0).")
  }

  wt_sequence <- toupper(wt_rows$nt_seq[[1]])
  wt_nucleotides <- strsplit(wt_sequence, "", fixed = TRUE)[[1]]

  variants <- df |>
    dplyr::filter(
      !is.na(Nham_nt),
      sequence_length == length(wt_nucleotides),
      !is.na(nt_seq),
      !is.na(relative_catalysis),
      !is.na(relative_leak),
      !is.na(net_catalysis)
    ) |>
    dplyr::mutate(
      .positions = purrr::map(nt_seq, \(sequence) {
        sequence_nucleotides <- strsplit(toupper(sequence), "", fixed = TRUE)[[1]]
        which(sequence_nucleotides != wt_nucleotides)
      })
    ) |>
    tidyr::unnest_longer(.positions, values_to = "position") |>
    dplyr::mutate(
      wt_nucleotide = wt_nucleotides[position],
      mutant_nucleotide = purrr::map_chr(
        nt_seq,
        \(sequence) strsplit(toupper(sequence), "", fixed = TRUE)[[1]][position]
      ),
      catalysis_change = relative_catalysis - 1,
      leak_change = relative_leak - 1,
      net_catalysis_change = net_catalysis - 1
    ) |>
    dplyr::filter(mutant_nucleotide %in% c("A", "C", "G", "T"))

  by_nucleotide <- variants |>
    dplyr::group_by(position, wt_nucleotide, mutant_nucleotide) |>
    dplyr::summarise(
      count = dplyr::n(),
      relative_catalysis_sum = sum(catalysis_change),
      relative_leak_sum = sum(leak_change),
      net_catalysis_sum = sum(net_catalysis_change),
      relative_catalysis_mean = mean(catalysis_change),
      relative_leak_mean = mean(leak_change),
      net_catalysis_mean = mean(net_catalysis_change),
      .groups = "drop"
    ) |>
    tidyr::complete(
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
    dplyr::mutate(wt_nucleotide = wt_nucleotides[position]) |>
    tidyr::pivot_wider(
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
      names_glue = "mutated_to_{mutant_nucleotide}_{.value}"
    )

  by_nucleotide |>
    dplyr::mutate(
      total_mutations = rowSums(
        dplyr::pick(dplyr::matches("^mutated_to_[ACGT]_count$"))
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      position,
      wt_nucleotide,
      total_mutations,
      dplyr::everything()
    )
}
