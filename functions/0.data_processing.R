#' Read and cleanup DiMSum output data
#' - Read data lazily (selective cols)
#' - cleanup: concise count column names
#' - calculate: frequency, enrichment and leak metrics
#' @param file_path path to the input file
#' @return a cleaned data frame

read_and_cleanup_data <- function(file_path) {

  # signatures: to detect cols = input, enrichment, leak 
  input_signature <- "Input"
  selection_signature <- "(neg25uM5fDUTheo|UraHisTheo|Enriched)"
  leak_signature <- "(NoTheo|UraHisVeh|Unenriched)"
  
  data <- read_tsv(file_path, lazy = TRUE) |>
    select("nt_seq", "Nham_nt", ends_with("_count")) |> # select only relevant columns

    # cleanup column names
    rename_with(
      .fn = ~ case_when(
        str_detect(.x, input_signature) ~ "Input_count", # detect signature and switch to concise name
        str_detect(.x, selection_signature) ~ "Selection_count",
        str_detect(.x, leak_signature) ~ "Leak_count",
        TRUE ~ .x  # Keep original name if no match
      ),
      .cols = ends_with("_count")
    ) |>

    # Add columns for calculations: sequence length, frequency of counts
    mutate(sequence_length = str_length(nt_seq)) |> # sequence length
    mutate(across(
      ends_with("_count"),
      ~ .x / sum(.x, na.rm = TRUE),
      .names = "{.col}_freq"
    )) |>
    rename_with(
      ~ str_replace(.x, "_count_freq$", "_freq") |> str_to_lower(),
      ends_with("_count_freq")
    ) %T>%
    
    { # Special section: not manipulating the data frame, but performing checks
      # error check for wild-type sequence (Nham_nt == 0); calculation continues in next step
      WT_count <- .$Nham_nt[.$Nham_nt == 0 & !is.na(.$Nham_nt)] %>% length() 
      
      # error check that breaks calculations
      if (WT_count == 1) {
        message(WT_count, " => A single wild-type sequence found (Nham_nt == 0).")
      } else {
        message("Cannot proceed with relative catalysis/leak calculations.")
        stop("Non unique/no wild-type sequences found: ", WT_count)
      }
    } |>

    # calculate enrichment and leak metrics
    mutate(
      catalysis_score = selection_freq / input_freq,
      leak_score = leak_freq / input_freq,
      relative_catalysis = catalysis_score / catalysis_score[Nham_nt == 0 & !is.na(Nham_nt)],
      relative_leak = leak_score / leak_score[Nham_nt == 0 & !is.na(Nham_nt)],
      net_catalysis = relative_catalysis / (1 + relative_leak)
  )
}