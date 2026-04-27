#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
})

parse_args <- function(args) {
  parsed <- list(
    input_tsv = NULL,
    out_prefix = NULL,
    dp_column = "mean_dp",
    title_prefix = "Sample DP"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1]] else NULL

    if (key == "--input-tsv") {
      parsed$input_tsv <- value
      i <- i + 2
    } else if (key == "--out-prefix") {
      parsed$out_prefix <- value
      i <- i + 2
    } else if (key == "--dp-column") {
      parsed$dp_column <- value
      i <- i + 2
    } else if (key == "--title-prefix") {
      parsed$title_prefix <- value
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  plot_sample_dp_distribution.R --input-tsv FILE --out-prefix PREFIX [--dp-column mean_dp] [--title-prefix STR]\n",
        sep = ""
      )
      quit(status = 0)
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  if (is.null(parsed$input_tsv) || is.null(parsed$out_prefix)) {
    stop("Both --input-tsv and --out-prefix are required.", call. = FALSE)
  }

  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

df <- read_tsv(args$input_tsv, show_col_types = FALSE)

if (!(args$dp_column %in% names(df))) {
  stop(sprintf("Column '%s' not found in %s", args$dp_column, args$input_tsv), call. = FALSE)
}

plot_df <- data.frame(dp_value = suppressWarnings(as.numeric(df[[args$dp_column]])))
plot_df <- plot_df[!is.na(plot_df$dp_value), , drop = FALSE]

if (nrow(plot_df) == 0) {
  stop("No valid DP values found for plotting.", call. = FALSE)
}

p_hist <- ggplot(plot_df, aes(x = dp_value)) +
  geom_histogram(bins = 60, fill = "grey70", color = "black") +
  labs(
    title = sprintf("%s distribution", args$title_prefix),
    x = "DP",
    y = "Number of samples"
  ) +
  theme_bw()

p_density <- ggplot(plot_df, aes(x = dp_value)) +
  geom_density(fill = "grey70", alpha = 0.5) +
  labs(
    title = sprintf("%s density", args$title_prefix),
    x = "DP",
    y = "Density"
  ) +
  theme_bw()

ggsave(sprintf("%s.histogram.pdf", args$out_prefix), p_hist, width = 8, height = 5)
ggsave(sprintf("%s.density.pdf", args$out_prefix), p_density, width = 8, height = 5)
