library(tidyverse)

#' Parses a dot-bracket notation string and returns a data frame with structural annotations.
#' Author: Claude.AI with Prashant K. (17/Nov/26)
#' @param dotbracket A string in dot-bracket notation representing RNA secondary structure.
#' @param pal A named vector of colors for different structural elements (stems, hairpins, interior loops, multiloops, external loops).
#' @return A data frame with columns: Name, arm, start, end, type, and colour, representing the structural elements of
parse_dotbracket <- function(dotbracket, pal = c(stem = "darkseagreen3", loop = "thistle3", tail = "grey80")) {
  ch <- strsplit(dotbracket, "")[[1]]; n <- length(ch) # split the dot-bracket string into characters and get its length
  
  # Initialize vectors to track partners, enclosing pairs, and parents for each dot-bracket char
  partner <- integer(n); enc <- integer(n); parent <- integer(n)
  stack <- integer() # stack to keep track of open parentheses

  # identify pairs and their relationships within the dotbracket string
  for (i in seq_len(n)) {
    top <- if (length(stack)) tail(stack, 1) else 0L
    if (ch[i] == "(") { parent[i] <- top; stack <- c(stack, i) }
    else if (ch[i] == ")") {
      j <- tail(stack, 1); stack <- head(stack, -1)
      partner[i] <- j; partner[j] <- i
    } else enc[i] <- top
  }

  # stems: runs of stacked pairs, one row per arm
  op  <- which(partner > seq_len(n)) # get indices of opening parentheses (5' end of stems)
  new <- c(TRUE, !(diff(op) == 1 & diff(partner[op]) == -1)) # identify new stems based on consecutive pairs
  
  # stems: create a data frame with start and end positions for each stem arm
  stems <- tibble(pos = op, pair = partner[op], id = cumsum(new)) |>
    summarise(`5p` = list(range(pos)), `3p` = list(range(pair)), .by = id) |>
    pivot_longer(c(`5p`, `3p`), names_to = "arm", values_to = "r") |>
    mutate(Name = paste0("stem", id), arm, start = map_int(r, 1), end = map_int(r, 2),
              type = "stem", .keep = "none")


  # loops: unpaired runs grouped by the enclosing pair

  # count the number of children for each enclosing pair (to classify loops: junctions vs closed loops)
  nchild <- table(factor(parent[op], levels = c(0, op))) # 
  runs <- rle(ch == ".")
  e <- cumsum(runs$lengths); s <- e - runs$lengths + 1
  
  # loops: create a data frame with start and end positions for each loop,
  # .. along with its classification (tail, junction, closed)
  loops <- tibble(start = s, end = e, dot = runs$values) |>
  filter(dot) |>
  mutate(enc   = enc[start],
          nc    = ifelse(enc == 0, NA, as.integer(nchild[as.character(enc)])),
          class = case_when(enc == 0 ~ "tail", nc >= 2 ~ "junction", TRUE ~ "closed")) |>
  mutate(grp = min(start), .by = c(class, enc)) |>
  mutate(Name = paste0(c(closed = "L", junction = "junction", tail = "tail")[class],
                        dense_rank(grp)), .by = class) |>
  mutate(Name, arm = NA_character_, start, end,
            type = ifelse(class == "tail", "tail", "loop"), .keep = "none")

  # Combine stems and loops, assign colors, and return the data frame
  bind_rows(stems, loops) |>
    arrange(start) |>
    mutate(colour = unname(pal[type])) |> 

    fill(arm, .direction = "downup") # fill arm column for loops based on the above stem
}

# Notes: 
# An error of this nature: Error in partner[i] <- j : replacement has length zero
# Suggests that the dot-bracket string is wrong; has unmatched parentheses. 
# Check the dot-bracket string for correctness.

# Testing 
# hhrzTheoAptaz_dotbracket <- ".(((((((((....((((((.......((((......)))).......(((((.(((((.......))))).))))).))))))....))))))))))."
# ann <- parse_dotbracket(hhrzTheoAptaz_dotbracket) |> print()

## Plot annotation ----

dark <- function(x) ifelse(is.na(x), "black", rgb(t(col2rgb(x) * 0.6), maxColorValue = 255))


#' Annotates a ggplot with RNA structural elements based on the provided annotation data frame.
#' Author: Claude.AI with Prashant K. (17/Nov/26)
#' @param ann_df A data frame containing structural annotations with columns: Name, arm, start, end, type, and colour.
#' @param pal A named vector of colors for different structural elements (stems, loops, tails).
#' @param vj A named vector specifying vertical justification for text annotations (auto and manual).
#' @param sz A named vector specifying text size for annotations (auto and manual).
#' @param a A named vector specifying alpha transparency for the rectangles (auto and manual).
#' @return A list of ggplot layers that can be added to a ggplot object to visualize the RNA structural elements.
annotate_sequence <- function(ann_df,
                         pal = c(stem = "darkseagreen3", loop = NA, tail = "grey80"),
                         vj  = c(auto = 3.5, manual = 1),   # label tier, relative to panel top
                         stagger = 1.4,  # vertical offset for auto labels
                         sz  = c(auto = 2.5, manual = 4),
                         a   = c(auto = 0.12, manual = 0.2),
                         face = c(auto = "plain", manual = "bold")) {
  ann_df <- ann_df |>
    filter(is.na(type) | type != "tail") |>                    # manual rows have no type
    mutate(colour = ifelse(source == "auto", unname(pal[type]), colour),
           txt = dark(colour), mid = (start + end) / 2,
           
           tier = ifelse(source == "auto", (cumsum(source == "auto") - 1) %% 3, 0), # prevent overlap by staggaring auto labels
           vj_ = ifelse(source == "auto", vj["auto"] + tier * stagger, vj["manual"]),
           sz_ = sz[source], a_ = a[source], face_ = face[source])

  # plot annotation rectangles and text labels
  # Note: I() plots the value as is, without scaling or transformation/legend
  list(
    geom_rect(data = ann_df, inherit.aes = FALSE,
              aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf,
                  fill = I(colour), alpha = I(a_))),
    geom_text(data = ann_df, inherit.aes = FALSE, y = Inf,
              aes(x = mid, label = Name, colour = I(txt), vjust = I(vj_),
                  size = I(sz_), fontface = I(face_))),
    coord_cartesian(clip = "off"),
    theme(plot.margin = margin(t = 25))
  )
}


# Plot contributions across the sequence -----
#' Plots the contribution of mutations across the ribozyme sequence, with optional annotations for structural elements.
#' Author: GH Copilot with Prashant K. (<17/Nov/26)
#' @param contribution_df A data frame containing contribution data with columns: position and the specified contribution column.
#' @param .col The name of the column in contribution_df to plot (default: "total_net_catalysis_mean").
#' @return A ggplot object visualizing the contribution across the sequence, with optional structural annotations.
plot_contribution_across_sequence <- 
  function(contribution_df, .col = "total_net_catalysis_mean",
           title = "Contribution assigned to mutations across sequence") {

  # plot total_net catalysis across positions
contribution_df |> 
  # filter(total_mutations > 0) |> 
  ggplot(aes(x = position, y = .data[[.col]])) +
  geom_point() + geom_line(alpha = 0.4) +

  labs(title = title) +
  theme_minimal()

}