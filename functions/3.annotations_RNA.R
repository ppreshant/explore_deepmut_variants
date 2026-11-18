library(tidyverse)

#' Parses a dot-bracket notation string and returns a data frame with structural annotations.
#' Author: Claude.AI with Prashant K. (17/Nov/26)
#' @param dotbracket A string in dot-bracket notation representing RNA secondary structure.
#' @param pal A named vector of colors for different structural elements (stems, hairpins, interior loops, multiloops, external loops).
#' @return A data frame with columns: Name, arm, start, end, type, and colour, representing the structural elements of
parse_dotbracket <- function(dotbracket, pal = c(stem = "#1b9e77", hairpin = "#d95f02",
                                         interior = "#7570b3", multiloop = "#e7298a",
                                         external = "#999999")) {
  ch <- strsplit(dotbracket, "")[[1]]; n <- length(ch)
  partner <- integer(n); enc <- integer(n); parent <- integer(n)
  stack <- integer()
  for (i in seq_len(n)) {
    top <- if (length(stack)) tail(stack, 1) else 0L
    if (ch[i] == "(") { parent[i] <- top; stack <- c(stack, i) }
    else if (ch[i] == ")") {
      j <- tail(stack, 1); stack <- head(stack, -1)
      partner[i] <- j; partner[j] <- i
    } else enc[i] <- top
  }

  # stems: consecutive stacked pairs
  op <- which(partner > seq_len(n))
  new <- c(TRUE, !(diff(op) == 1 & diff(partner[op]) == -1))
  stems <- tibble(pos = op, pair = partner[op], id = cumsum(new)) |>
    group_by(id) |>
    summarise(`5p` = list(range(pos)), `3p` = list(range(pair)), .groups = "drop") |>
    pivot_longer(c(`5p`, `3p`), names_to = "arm", values_to = "r") |>
    transmute(Name = paste0("stem", id), arm, start = map_int(r, 1), end = map_int(r, 2),
              type = "stem")

  # loops: runs of unpaired bases, typed by the enclosing pair's number of child helices
  nchild <- table(factor(parent[op], levels = c(0, op)))
  runs <- rle(ch == ".")
  e <- cumsum(runs$lengths); s <- e - runs$lengths + 1
  loops <- tibble(start = s, end = e, dot = runs$values) |>
    filter(dot) |>
    mutate(enc = enc[start],
           nc = ifelse(enc == 0, NA, as.integer(nchild[as.character(enc)])),
           type = case_when(enc == 0 ~ "external", nc == 0 ~ "hairpin",
                            nc == 1 ~ "interior", TRUE ~ "multiloop")) |>
    group_by(type) |> mutate(Name = paste0(type, row_number())) |> ungroup() |>
    mutate(arm = NA_character_) |>
    select(Name, arm, start, end, type)

  bind_rows(stems, loops) |>
    arrange(start) |>
    mutate(colour = unname(pal[type]))
}

# Notes: 
# An error of this nature: Error in partner[i] <- j : replacement has length zero
# Suggests that the dot-bracket string is wrong; has unmatched parentheses. 
# Check the dot-bracket string for correctness.

# Testing 
# hhrzTheoAptaz_dotbracket <- ".(((((((((....((((((.......((((......)))).......(((((.(((((.......))))).))))).))))))....))))))))))."
# ann <- parse_dotbracket(hhrzTheoAptaz_dotbracket) |> print()