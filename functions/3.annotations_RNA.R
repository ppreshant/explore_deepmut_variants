library(tidyverse)

#' Parses a dot-bracket notation string and returns a data frame with structural annotations.
#' Author: Claude.AI with Prashant K. (17/Nov/26)
#' @param dotbracket A string in dot-bracket notation representing RNA secondary structure.
#' @param pal A named vector of colors for different structural elements (stems, hairpins, interior loops, multiloops, external loops).
#' @return A data frame with columns: Name, arm, start, end, type, and colour, representing the structural elements of
parse_dotbracket <- function(db, pal = c(stem = "#1b9e77", loop = "#d95f02", tail = "#999999")) {
  ch <- strsplit(db, "")[[1]]; n <- length(ch)
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

  # stems: runs of stacked pairs, one row per arm
  op  <- which(partner > seq_len(n))
  new <- c(TRUE, !(diff(op) == 1 & diff(partner[op]) == -1))
  stems <- tibble(pos = op, pair = partner[op], id = cumsum(new)) |>
        summarise(`5p` = list(range(pos)), `3p` = list(range(pair)), .by = id) |>
    pivot_longer(c(`5p`, `3p`), names_to = "arm", values_to = "r") |>
    transmute(Name = paste0("stem", id), arm, start = map_int(r, 1), end = map_int(r, 2),
              type = "stem")

  # loops: unpaired runs grouped by the enclosing pair
  nchild <- table(factor(parent[op], levels = c(0, op)))
  runs <- rle(ch == ".")
  e <- cumsum(runs$lengths); s <- e - runs$lengths + 1

  # loops: unpaired runs grouped by the enclosing pair
  loops <- tibble(start = s, end = e, dot = runs$values) |>
    filter(dot) |>
    mutate(enc   = enc[start],
           nc    = ifelse(enc == 0, NA, as.integer(nchild[as.character(enc)])),
           class = case_when(enc == 0 ~ "tail", nc >= 2 ~ "junction", TRUE ~ "closed")) |>
    mutate(grp = min(start), .by = c(class, enc)) |>
mutate(Name = paste0(c(closed = "L", junction = "junction", tail = "tail")[class],
                        dense_rank(grp)), .by = class) |>
  transmute(Name, arm = NA_character_, start, end,
            type = ifelse(class == "tail", "tail", "loop"))

  # Combine stems and loops, assign colors, and return the data frame
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