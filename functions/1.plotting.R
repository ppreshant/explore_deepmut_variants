# Plotting functions ---------

## annotation generating function -----

#' Generates a function to create text annotations for a plot.
#' @param max_x: the maximum x value in the plot (used to position the annotation)
#' @param max_y: the maximum y value in the plot (used to position the annotation)
#' @return: A function that takes x_scale, label_text, and y_offset as arguments and returns a ggplot2 annotation layer.

annotation_generator <- function(max_x_val, max_y_val) { 

  function(x_scale = 0.6, label_text, y_offset) {
  # work in progress ;
  # make a text annotation to add to the plot; based on offset values for x and y
    annotate(geom = "text",
      x = max_x_val * x_scale, 
      y = 0 + y_offset * max_y_val/10, 
      label = label_text, 
      hjust = 0, vjust = 1.5, size = 3)
  }
}


## plot counts vs nham -----
#' Plot the distribution of Nham_nt counts in a data frame.
#' @param df A data frame containing a column named Nham_nt.
#' @param print_stats Whether to print summary statistics.
#' @return A ggplot object showing the distribution of Nham_nt counts.

nham_any_count_distribution <-  function(
    counts,
    x_identifier = "variants",
    print_stats = FALSE) {

  # TODO: format the axis labels and annotations for pretty reading (comman in numbers etc.)

  # calculate summaries for plotting
  snp_count <- counts |> filter(Nham_nt > 0) |> pull(count) |> sum()
  indel_count <- counts |> filter(is.na(Nham_nt)) |> pull(count) |> sum()
  total_count <- counts |> pull(count) |> sum()

  # Calculate dynamic positions
  max_count <- max(counts$count, na.rm = FALSE)
  num_nham <- nrow(counts)
  add_annotation <- annotation_generator(max_x_val = max_count, max_y_val = num_nham) # create a function to add annotations to the plot
 
  if(print_stats) {
    # print stats
    cat("snp_count =", snp_count, "\nindel_count = ", indel_count, "\n")
    cat("total = ", total_count, "variants\n")
  }

  # Plot variant counts
  plt <- ggplot(counts, aes(count, fct_rev(as_factor(Nham_nt)))) + 
    geom_bar(stat = "identity", width = 1) + 

    # annotations
    labs(x = str_c("Number of ", x_identifier), y = "Number of mismatches") + # axis labels
    geom_text(aes(label = format_numeric_labels(count)), hjust = -0.1, size = 3) + # count of each bar
    # summary variant counts: snps, indels, total
    add_annotation(label_text = str_c("snps: ", format_numeric_labels(snp_count)), y_offset = 2) +
    add_annotation(label_text = str_c("indels: ", format_numeric_labels(indel_count)), y_offset = 1.5) +
    add_annotation(label_text = str_c("total: ", format_numeric_labels(total_count)), y_offset = 1) +

    # aesthetics
    theme_bw() + 
    scale_x_continuous(expand = expansion(mult = c(0.1, .25))) # expand x-axis to make room for text annotations

  # print(plt)
  
}

## Read distribution -----
read_distribution <- function(df, .col = "Input_count") {

  ggplot(df, aes(y = .data[[.col]])) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    scale_y_log10() + 
    scale_x_continuous(name = "Count (thousands)", 
      labels = scales::label_number(scale = 1e-3)) +

    # annotations
    # stat_bin( # add counts on top of bars
    #   binwidth = 5,
    #   geom = "text",
    #   aes(label = after_stat(count)), 
    #   vjust = -0.4) +
    labs(y = .col) +
    theme_minimal()

}

### wrapper to patch this across columns 

patch_read_dist_multi_columns <- 
  function(
    df, 
    .cols = c("Input_count", "Selection_count", "Leak_count"),
    title = "Distribution of counts by column"
) {

  # Create a list of plots for each column
  plots <- 
    lapply(.cols, read_distribution, df = df) |>
      patchwork::wrap_plots(ncol = length(.cols), axes = "collect")
    
  plots + patchwork::plot_annotation(title = title)
}

# Summarize and send to plot -----

## plot variant counts vs nham -----
#' Plot the distribution of Nham_nt counts in a data frame.
#' @param df A data frame containing a column named Nham_nt.
#' @param print_stats Whether to print summary statistics.
#' @return A ggplot object showing the distribution of Nham_nt counts.
 
nham_variant_count_distribution <-  function(df, print_stats = FALSE) {

  counts <- summarise(df, count = n(), .by = Nham_nt) # summarize to get variant counts at each Nham_nt value
  nham_any_count_distribution(counts) # call the generic function to plot counts vs nham
}

## plot read counts vs nham -----
#' Plot the distribution of reads by Nham_nt.
#' @param df A data frame containing a column named Nham_nt.
#' @param print_stats Whether to print summary statistics.
#' @return A ggplot object showing the distribution of Nham_nt counts.

nham_read_count_distribution <- function(df, print_stats = FALSE) {
  
  counts <- summarise(df, count = sum(Input_count), .by = Nham_nt) # summarize to get variant counts at each Nham_nt value
  nham_any_count_distribution(counts, x_identifier = "reads") # call the generic function to plot read counts vs nham
}

variant_count_by_threshold <- function(df, filtered_columns, 
                                      expected_count = expected_variants) {

  # plot the variant counts vs expected variants from count_tracker
  ggplot(df, aes(x = threshold, y = variants, label = threshold)) +
    geom_line(alpha = 0.7) +
    geom_point() + 
    ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf, alpha = 0.7) + # add labels to points

    # expected variants line
    geom_hline(yintercept = expected_count, linetype = "dashed") +
    annotate("text", 
      x = max(df$threshold) * 0.7, 
      y = expected_count + max(df$variants)/10, 
      label = str_c("Expected: ", format(expected_count, big.mark = ","), " variants")) +
    
    # annotations
    labs(title = "Variant Counts with thresholds",
        subtitle = str_c("Filtered columns: ", str_c(filtered_columns, collapse = ", ")),
        x = "Count Threshold",
        y = "Number of Variants") +
    theme_minimal()

}

## Plot contributions across the sequence -----
plot_contribution_across_sequence <- 
  function(contribution_df, .col = "total_net_catalysis_mean",
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


# formatting numbers for plot labels ----
#' Format numbers for plot labels: use scientific notation for large numbers, and commas for thousands.
#' @param n A numeric vector to format.
#' @return A character vector with formatted numbers.
format_numeric_labels <- function(n) {
  # Format numbers for plot labels: use scientific notation for large numbers, and commas for thousands
  case_when(
    n >= 1e6 ~ lab_sci(n),  # ex: 2,000,000 -> 2e6
    n >= 1e3 ~ str_c(format(round(n / 1000), big.mark = ","), "k"),  # ex: 2,000 -> 2K
    n >= 1 ~ format(n, scientific = FALSE),  # format for numbers >= 1
    TRUE ~ format(n, digits = 2)  # default case for numbers < 1
  )
}

## format numbers for plot labels in scientific notation ----
#' Format numbers for plot labels in scientific notation, removing unnecessary leading zeros.
#' @param x A numeric vector to format.
#' @return A character vector with formatted numbers.
lab_sci <- function(x) {
  sub("e\\+?(-?)0*", "e\\1",
      scales::label_scientific(digits = 1)(x))
}