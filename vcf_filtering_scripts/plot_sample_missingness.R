#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
})

parse_args <- function(args) {
  parsed <- list(
    input_tsv = NULL,
    out_prefix = NULL,
    cutoffs_pct = c(5, 10, 20, 40),
    title_prefix = "Sample missingness"
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
    } else if (key == "--cutoffs-pct") {
      parsed$cutoffs_pct <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]))
      i <- i + 2
    } else if (key == "--title-prefix") {
      parsed$title_prefix <- value
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  plot_sample_missingness.R --input-tsv FILE --out-prefix PREFIX [--cutoffs-pct 5,10,20,40] [--title-prefix STR]\n",
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

  parsed$cutoffs_pct <- parsed$cutoffs_pct[!is.na(parsed$cutoffs_pct)]
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

df <- read_tsv(args$input_tsv, show_col_types = FALSE)

required_cols <- c("sample", "missing_rate")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(
    sprintf("Missing required column(s) in %s: %s", args$input_tsv, paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

df$missing_rate <- suppressWarnings(as.numeric(df$missing_rate))
df <- df[!is.na(df$missing_rate), , drop = FALSE]

if (nrow(df) == 0) {
  stop("No valid missing_rate values found for plotting.", call. = FALSE)
}

df$missing_pct <- df$missing_rate * 100
cutoffs <- seq(0, 1, by = 0.01)
cut_df <- data.frame(
  cutoff = cutoffs,
  n_kept = vapply(cutoffs, function(c) sum(df$missing_rate <= c, na.rm = TRUE), numeric(1))
)
cut_df$cutoff_pct <- cut_df$cutoff * 100

p_hist <- ggplot(df, aes(x = missing_pct)) +
  geom_histogram(bins = 60, fill = "grey70", color = "black") +
  geom_vline(xintercept = args$cutoffs_pct, linetype = 2, color = "red") +
  labs(
    title = sprintf("%s distribution", args$title_prefix),
    x = "Missingness (%)",
    y = "Number of samples"
  ) +
  theme_bw()

p_density <- ggplot(df, aes(x = missing_pct)) +
  geom_density(fill = "grey70", alpha = 0.5) +
  geom_vline(xintercept = args$cutoffs_pct, linetype = 2, color = "red") +
  labs(
    title = sprintf("%s density", args$title_prefix),
    x = "Missingness (%)",
    y = "Density"
  ) +
  theme_bw()

p_cutoff <- ggplot(cut_df, aes(x = cutoff_pct, y = n_kept)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_vline(xintercept = args$cutoffs_pct, linetype = 2, color = "red") +
  labs(
    title = sprintf("%s cutoff vs samples kept", args$title_prefix),
    x = "Missingness cutoff (%)",
    y = "Number of samples kept"
  ) +
  theme_bw()

ggsave(sprintf("%s.histogram.pdf", args$out_prefix), p_hist, width = 8, height = 5)
ggsave(sprintf("%s.density.pdf", args$out_prefix), p_density, width = 8, height = 5)
ggsave(sprintf("%s.cutoff_vs_samples_kept.pdf", args$out_prefix), p_cutoff, width = 8, height = 5)
