library(tidyverse)

#' Parses a dot-bracket notation string and returns a data frame with structural annotations.
#' Author: Claude.AI with Prashant K. (17/Nov/26)
#' @param dotbracket A string in dot-bracket notation representing RNA secondary structure.
#' @param pal A named vector of colors for different structural elements (stems, hairpins, interior loops, multiloops, external loops).
#' @return A data frame with columns: Name, arm, start, end, type, and colour, representing the structural elements of
parse_dotbracket <- function(db, pal = c(stem = "darkseagreen3", loop = "thistle3", tail = "grey80")) {
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


# Plot contributions across the sequence -----
#' Plots the contribution of mutations across the ribozyme sequence, with optional annotations for structural elements.
#' Author: GH Copilot with Prashant K. (<17/Nov/26)
#' @param contribution_df A data frame containing contribution data with columns: position and the specified contribution column.
#' @param .col The name of the column in contribution_df to plot (default: "total_net_catalysis_mean").
#' @param annotation_df An optional data frame containing structural annotations with columns: Name, arm, start, end, type, and colour.
#' @return A ggplot object visualizing the contribution across the sequence, with optional structural annotations.
plot_contribution_across_sequence <- 
  function(contribution_df, .col = "total_net_catalysis_mean",
           annotation_df = NULL,
           title = "Contribution assigned to mutations across sequence") {

  # plot total_net catalysis across positions
contribution_df |> 
  # filter(total_mutations > 0) |> 
  ggplot(aes(x = position, y = .data[[.col]])) +
  geom_point() + geom_line(alpha = 0.4) +

  # annotations for positions of interest
  geom_hline(yintercept = 0.5, color = "red") + 
  annotate("rect", xmin = 21, xmax = 27, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "orange") + # ribozyme core
  annotate("text", x = 24, y = 2.5, label = "Core", vjust = -1, color = "darkorange3") +

  # communication modules
  annotate("rect", xmin = 49, xmax = 53, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
  annotate("text", x = 51, y = 2.5, label = "CM1", vjust = -1, color = "blue") +

  annotate("rect", xmin = 73, xmax = 77, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
  annotate("text", x = 75, y = 2.5, label = "CM2", vjust = -1, color = "blue") +
  
  # aptamer
  annotate("text", x = 63, y = 2.5, label = "Aptamer", vjust = -1) +

  labs(title = title) +
  theme_minimal()

}