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
    geom_text(aes(label = count), hjust = -0.1, size = 3) + # count of each bar
    # summary variant counts: snps, indels, total
    add_annotation(label_text = str_c("snps: ", snp_count), y_offset = 2) +
    add_annotation(label_text = str_c("indels: ", indel_count), y_offset = 1.5) +
    add_annotation(label_text = str_c("total: ", total_count), y_offset = 1) +

    # aesthetics
    theme_bw() + 
    scale_x_continuous(expand = expansion(mult = c(0.1, .25))) # expand x-axis to make room for text annotations
  print(plt)
}

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