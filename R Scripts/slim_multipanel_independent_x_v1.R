#!/usr/bin/env Rscript
# slim_multipanel_generator.R
# -----------------------------------------------------------------------------
# Combines 2 or 4 SLiM simulation outputs into a single multi-panel figure,
# styled to match the Rupe (JoMB) manuscript figures (v9 styling).
#
# Input files can be EITHER:
#   - All trajectory logs (CSV header begins with "cycle,Population_size")
#       -> Multi-panel fitness trajectory comparison (like Fig. 14)
#   - All full-output dumps (first line begins with "// Initial random seed:")
#       -> Multi-panel allele frequency histogram comparison
# Mixing types within one figure is rejected.
#
# Usage:
#   Rscript slim_multipanel_generator.R <file1> <file2> [file3 file4] [--format png|pdf|...]
#
# Panel labels (a, b, c, d) are added automatically in the top-left of each panel
# based on input order. A single shared legend is placed at the top of the figure.
# All panels use the same y-axis range so they are directly comparable.
#
# Layouts:
#   2 files -> side-by-side (1 row, 2 columns), 14 x 5.5 in
#   4 files -> 2x2 grid, 14 x 10 in
#
# Output: written next to the FIRST input file as
#   <basename_of_first>-multipanel-<timestamp>.<ext>
#
# Per-panel descriptions belong in the manuscript caption (e.g., "(a) H&C baseline;
# (b) B1 Run 1; ..."), not on the figure itself.
# -----------------------------------------------------------------------------

# ---------- dependency loading ------------------------------------------------
require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Missing R package '%s'. Install with: install.packages(\"%s\")",
                 pkg, pkg), call. = FALSE)
  }
}
require_pkg("optparse")
require_pkg("ggplot2")
require_pkg("scales")
require_pkg("patchwork")   # multi-panel figure composition

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ---------- JoMB style constants (v9) -----------------------------------------
JOMB_NAVY      <- "#08306B"
JOMB_CRIMSON   <- "#B1182C"
JOMB_FONT      <- "Arial"
JOMB_BASE_SIZE <- 15
JOMB_LINE_W    <- 0.7   # trajectory line width

theme_jomb <- function(base_size = JOMB_BASE_SIZE,
                       base_family = JOMB_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title.x = element_text(face = "bold", size = base_size,
                                  margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = base_size,
                                  margin = margin(r = 8)),
      axis.text.x  = element_text(size = base_size - 1, color = "black"),
      axis.text.y  = element_text(size = base_size - 1, color = "black"),
      axis.line    = element_line(color = "black", size = 0.6),
      axis.ticks   = element_line(color = "black", size = 0.6),
      axis.ticks.length = unit(4, "pt"),

      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),

      # In multi-panel mode each panel's legend is suppressed (we add a single
      # shared legend at the figure level via patchwork). But we keep these
      # settings so individual panels can still render correctly if needed.
      legend.position      = c(0.03, 0.97),
      legend.justification = c("left", "top"),
      legend.background    = element_blank(),
      legend.key           = element_rect(fill = "white", color = NA),
      legend.title         = element_blank(),
      legend.text          = element_text(size = base_size - 1),
      legend.spacing.y     = unit(2, "pt"),
      legend.margin        = margin(2, 4, 2, 4),

      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin      = margin(10, 14, 6, 10),
      plot.title       = element_text(hjust = 0.5,
                                      size = base_size + 1,
                                      margin = margin(b = 8)),

      # Panel tag (the "a", "b", "c", "d" label) placement and styling
      plot.tag          = element_text(face = "bold", size = base_size + 4,
                                       family = JOMB_FONT),
      plot.tag.position = c(0.02, 0.97)
    )
}

# ---------- CLI parsing -------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--format"),
              type = "character",
              default = "png",
              help = "Output image format: png (default), pdf, jpeg, tiff, svg",
              metavar = "FORMAT")
)
parser <- OptionParser(
  usage = paste("Usage: %prog <file1> <file2> [<file3> <file4>] ",
                "[--format png|pdf|jpeg|tiff|svg]"),
  option_list = option_list
)
# `positional_arguments = c(2, 4)` allows either 2 or 4 positional args, but
# optparse's range checking is exclusive on the upper bound in some versions,
# so we accept up to 10 here and validate manually below.
parsed <- tryCatch(
  parse_args(parser, positional_arguments = c(2, 10)),
  error = function(e) {
    message("Argument error: ", conditionMessage(e))
    print_help(parser)
    quit(status = 2)
  }
)
input_paths <- parsed$args
output_format <- tolower(parsed$options$format)

# Strict validation: exactly 2 or 4 files allowed
if (!(length(input_paths) %in% c(2, 4))) {
  stop(sprintf(
    "This script requires exactly 2 or 4 input files; got %d.",
    length(input_paths)), call. = FALSE)
}

allowed_formats <- c("png", "pdf", "jpeg", "jpg", "tiff", "svg")
if (!(output_format %in% allowed_formats)) {
  stop(sprintf("Unsupported --format '%s'. Allowed: %s",
               output_format, paste(allowed_formats, collapse = ", ")),
       call. = FALSE)
}

# ---------- helpers (shared) --------------------------------------------------
detect_type <- function(path) {
  first_line <- tryCatch(
    readLines(path, n = 1, warn = FALSE),
    error = function(e) character(0)
  )
  if (length(first_line) == 0) return("unknown")
  if (startsWith(first_line, "cycle,Population_size")) return("trajectory")
  if (startsWith(first_line, "// Initial random seed:")) return("fullout")
  return("unknown")
}

require_columns <- function(df, needed) {
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Trajectory file is missing required column(s): %s.\nColumns found: %s",
      paste(missing, collapse = ", "),
      paste(colnames(df), collapse = ", ")
    ), call. = FALSE)
  }
}

make_output_path <- function(first_input_path, ext, n_panels) {
  dir <- dirname(first_input_path)
  stem <- tools::file_path_sans_ext(basename(first_input_path))
  ts <- format(Sys.time(), "%Y_%m_%d-%H.%M.%S")
  file.path(dir,
            sprintf("%s-%dpanel-%s.%s", stem, n_panels, ts, ext))
}

# ---------- per-panel builders ------------------------------------------------
# These return ggplot objects WITHOUT a legend (legend stripped here and
# re-added at the figure level via patchwork's plot_layout(guides = "collect")).

build_trajectory_panel <- function(input_path, y_max, x_max) {
  df <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(df, c("cycle", "Mean_fitness", "Intrinsic_fitness"))
  gen_kx <- df[["cycle"]] / 1000
  long_df <- rbind(
    data.frame(gen_kx = gen_kx, fitness = df[["Mean_fitness"]],
               series = "Mean fitness (density-scaled)",
               stringsAsFactors = FALSE),
    data.frame(gen_kx = gen_kx, fitness = df[["Intrinsic_fitness"]],
               series = "Intrinsic fitness (genomic)",
               stringsAsFactors = FALSE)
  )
  long_df$series <- factor(
    long_df$series,
    levels = c("Mean fitness (density-scaled)", "Intrinsic fitness (genomic)")
  )
  color_map <- c(
    "Mean fitness (density-scaled)" = JOMB_NAVY,
    "Intrinsic fitness (genomic)"   = JOMB_CRIMSON
  )
  ggplot(long_df, aes(x = gen_kx, y = fitness, color = series)) +
    geom_line(linewidth = JOMB_LINE_W) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_x_continuous(
      name = expression(bold(paste("Generation (", "\u00D7", "1,000)"))),
      limits = c(0, x_max),
      expand = expansion(mult = c(0.005, 0.02))
    ) +
    scale_y_continuous(
      name = expression(bold("Fitness")),
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0))
    ) +
    theme_jomb()
}

# Single-pass parser for SLiM full-output mutation + genome sections.
# Returns: list(del_freqs, ben_freqs) — vectors of allele frequencies
# for deleterious (m2) and beneficial (m3) mutations.
parse_fullout_frequencies <- function(input_path) {
  lines <- readLines(input_path, warn = FALSE)
  section <- "header"
  mut_ids_acc   <- character(0)
  mut_types_acc <- character(0)
  genome_count <- 0L
  genome_ids_chunks <- list()
  for (raw in lines) {
    if (raw == "Mutations:")        { section <- "mutations";   next }
    if (raw == "Individuals:")      { section <- "individuals"; next }
    if (raw == "Genomes:")          { section <- "genomes";     next }
    if (startsWith(raw, "Populations:")) { section <- "populations"; next }
    if (section %in% c("header", "populations", "individuals")) next
    if (!nzchar(raw)) next
    if (section == "mutations") {
      fields <- strsplit(raw, "\\s+")[[1]]
      if (length(fields) >= 3) {
        mut_ids_acc[length(mut_ids_acc) + 1L]     <- fields[1]
        mut_types_acc[length(mut_types_acc) + 1L] <- fields[3]
      }
      next
    }
    if (section == "genomes") {
      if (!startsWith(raw, "p1:")) next
      genome_count <- genome_count + 1L
      fields <- strsplit(raw, "\\s+")[[1]]
      if (length(fields) > 2) {
        ids <- suppressWarnings(as.integer(fields[-c(1, 2)]))
        ids <- ids[!is.na(ids)]
        genome_ids_chunks[[length(genome_ids_chunks) + 1L]] <- ids
      }
      next
    }
  }
  if (genome_count == 0L) {
    stop(sprintf("No genome lines found in '%s' - is this a full-output file?",
                 input_path), call. = FALSE)
  }
  if (length(mut_ids_acc) == 0L) {
    stop(sprintf("No mutations parsed from '%s'.", input_path), call. = FALSE)
  }
  mut_type_map <- mut_types_acc
  names(mut_type_map) <- mut_ids_acc
  all_ids <- unlist(genome_ids_chunks, use.names = FALSE)
  counts <- table(all_ids)
  freqs <- as.numeric(counts) / genome_count
  freq_ids <- names(counts)
  types_for_freq <- mut_type_map[freq_ids]
  list(
    del = freqs[!is.na(types_for_freq) & types_for_freq == "m2"],
    ben = freqs[!is.na(types_for_freq) & types_for_freq == "m3"]
  )
}

build_histogram_panel <- function(input_path, y_max) {
  freqs <- parse_fullout_frequencies(input_path)
  del_freqs <- freqs$del
  ben_freqs <- freqs$ben

  FIXED_THRESHOLD <- 0.999
  BIN_WIDTH       <- 0.05
  del_fixed <- sum(del_freqs >= FIXED_THRESHOLD)
  ben_fixed <- sum(ben_freqs >= FIXED_THRESHOLD)
  del_nonfixed <- del_freqs[del_freqs < FIXED_THRESHOLD]
  ben_nonfixed <- ben_freqs[ben_freqs < FIXED_THRESHOLD]
  edges <- seq(0, 1, by = BIN_WIDTH)
  bin_centers <- edges[-length(edges)] + BIN_WIDTH / 2
  bin_idx_del <- cut(del_nonfixed, breaks = edges,
                     include.lowest = TRUE, right = FALSE, labels = FALSE)
  bin_idx_ben <- cut(ben_nonfixed, breaks = edges,
                     include.lowest = TRUE, right = FALSE, labels = FALSE)
  n_bins <- length(bin_centers)
  del_counts <- tabulate(bin_idx_del, nbins = n_bins)
  ben_counts <- tabulate(bin_idx_ben, nbins = n_bins)
  zero_to_na <- function(v) ifelse(v == 0L, NA_integer_, v)
  del_label <- sprintf("Deleterious (fixed: %d)", del_fixed)
  ben_label <- sprintf("Beneficial (fixed: %d)", ben_fixed)
  FIXED_X <- 1.0 + 2 * BIN_WIDTH
  plot_df <- rbind(
    data.frame(bin_x = bin_centers, count = zero_to_na(del_counts),
               category = del_label, stringsAsFactors = FALSE),
    data.frame(bin_x = bin_centers, count = zero_to_na(ben_counts),
               category = ben_label, stringsAsFactors = FALSE),
    data.frame(bin_x = FIXED_X,
               count = ifelse(del_fixed == 0L, NA_integer_, del_fixed),
               category = del_label, stringsAsFactors = FALSE),
    data.frame(bin_x = FIXED_X,
               count = ifelse(ben_fixed == 0L, NA_integer_, ben_fixed),
               category = ben_label, stringsAsFactors = FALSE)
  )
  plot_df$category <- factor(plot_df$category,
                             levels = c(del_label, ben_label))
  fill_map <- setNames(c(JOMB_CRIMSON, JOMB_NAVY), c(del_label, ben_label))
  bar_width   <- BIN_WIDTH * 0.4
  dodge_width <- BIN_WIDTH * 0.85
  tick_positions <- c(seq(0, 1, by = 0.1), FIXED_X)
  tick_labels    <- c(sprintf("%.1f", seq(0, 1, by = 0.1)), "Fixed")
  ggplot(plot_df, aes(x = bin_x, y = count, fill = category)) +
    geom_col(position = position_dodge(width = dodge_width),
             width = bar_width, na.rm = TRUE, color = NA) +
    geom_vline(xintercept = 1.0 + BIN_WIDTH / 2,
               linetype = "dotted", color = "black", linewidth = 0.4) +
    scale_fill_manual(values = fill_map, name = NULL) +
    scale_y_log10(
      name = expression(bold("Number of mutations (log scale)")),
      labels = trans_format("log10", math_format(10^.x)),
      limits = c(NA, y_max),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(
      name = expression(bold("Allele frequency")),
      breaks = tick_positions,
      labels = tick_labels,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    theme_jomb()
}

# ---------- main --------------------------------------------------------------
main <- function(input_paths, ext) {
  # Validate every file exists and is a .txt
  for (p in input_paths) {
    if (!file.exists(p))
      stop(sprintf("Path does not exist: %s", p), call. = FALSE)
    if (dir.exists(p))
      stop(sprintf("Folders are not accepted: %s", p), call. = FALSE)
    if (tolower(tools::file_ext(p)) != "txt")
      stop(sprintf("Input must have .txt extension: %s", p), call. = FALSE)
  }

  # Detect type for each file; reject mixed
  kinds <- sapply(input_paths, detect_type)
  if (any(kinds == "unknown")) {
    bad <- input_paths[kinds == "unknown"]
    stop(sprintf("Unrecognized SLiM file format(s):\n  %s",
                 paste(bad, collapse = "\n  ")), call. = FALSE)
  }
  if (length(unique(kinds)) != 1L) {
    stop("Mixed file types: all input files must be either all trajectory ",
         "logs or all full-output dumps, not a mix.", call. = FALSE)
  }
  kind <- kinds[[1]]
  n_panels <- length(input_paths)

  # ---- Build panels with shared axis ranges so panels are comparable ----
  if (kind == "trajectory") {
    # First pass: find each panel's OWN x_max (so panels with shorter runs don't
    # get padded with empty space out to the longest run). Y-axis is shared
    # across panels so the fitness values remain directly comparable.
    panel_x_max <- numeric(length(input_paths))
    all_y_max <- 1.1
    for (i in seq_along(input_paths)) {
      df <- read.csv(input_paths[[i]], check.names = FALSE, stringsAsFactors = FALSE)
      require_columns(df, c("cycle", "Mean_fitness", "Intrinsic_fitness"))
      panel_x_max[i] <- max(df[["cycle"]] / 1000, na.rm = TRUE)
      all_y_max <- max(all_y_max,
                       max(df[["Mean_fitness"]], na.rm = TRUE) * 1.05,
                       max(df[["Intrinsic_fitness"]], na.rm = TRUE) * 1.05)
    }
    # Build each panel with its own x_max but shared y_max
    panels <- lapply(seq_along(input_paths), function(i) {
      build_trajectory_panel(input_paths[[i]],
                             y_max = all_y_max,
                             x_max = panel_x_max[i])
    })
  } else {
    # Histograms: find the global max count across all panels so y-axis matches.
    all_y_max <- 1
    for (p in input_paths) {
      freqs <- parse_fullout_frequencies(p)
      # Quick max across the largest bin counts
      bw <- 0.05; thresh <- 0.999
      del_nf <- freqs$del[freqs$del < thresh]
      ben_nf <- freqs$ben[freqs$ben < thresh]
      n_bins <- length(seq(0, 1, by = bw)) - 1
      del_c <- if (length(del_nf) > 0)
        max(tabulate(cut(del_nf, breaks = seq(0, 1, by = bw),
                         include.lowest = TRUE, right = FALSE,
                         labels = FALSE), nbins = n_bins))
        else 0
      ben_c <- if (length(ben_nf) > 0)
        max(tabulate(cut(ben_nf, breaks = seq(0, 1, by = bw),
                         include.lowest = TRUE, right = FALSE,
                         labels = FALSE), nbins = n_bins))
        else 0
      df_fixed <- sum(freqs$del >= thresh)
      bf_fixed <- sum(freqs$ben >= thresh)
      panel_max <- max(del_c, ben_c, df_fixed, bf_fixed)
      all_y_max <- max(all_y_max, panel_max)
    }
    # Pad upward for log-scale headroom
    all_y_max <- all_y_max * 1.5
    panels <- lapply(input_paths, build_histogram_panel, y_max = all_y_max)
  }

  # ---- Strip legend from individual panels; patchwork will add a shared one ----
  panels <- lapply(panels, function(p) p + theme(legend.position = "none"))

  # ---- Add panel tags (a, b, c, d) ----
  tag_letters <- letters[seq_len(n_panels)]   # "a", "b", "c", "d"
  for (i in seq_len(n_panels)) {
    panels[[i]] <- panels[[i]] + labs(tag = tag_letters[i])
  }

  # ---- Compose with patchwork ----
  # For the shared legend: build ONE panel that retains its legend, then
  # extract that legend and place it at top of the combined figure.
  # Easier approach: use plot_layout(guides = "collect") with a SINGLE panel
  # that keeps its legend layer, then suppress it on the others. patchwork
  # handles guide-collection across all panels.
  #
  # We'll re-add the legend mapping to panel 1 only.
  legend_panel <- if (kind == "trajectory") {
    panels[[1]] + theme(legend.position = "right")
  } else {
    panels[[1]] + theme(legend.position = "right")
  }
  panels[[1]] <- legend_panel

  if (n_panels == 2) {
    combined <- panels[[1]] + panels[[2]] +
      plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "top",
            legend.justification = "left",
            legend.direction = "vertical",
            legend.box.margin = margin(0, 0, 6, 0))
    fig_w <- 14; fig_h <- 5.5
  } else {
    combined <- (panels[[1]] + panels[[2]]) /
                (panels[[3]] + panels[[4]]) +
      plot_layout(guides = "collect") &
      theme(legend.position = "top",
            legend.justification = "left",
            legend.direction = "vertical",
            legend.box.margin = margin(0, 0, 6, 0))
    fig_w <- 14; fig_h <- 10
  }

  out_path <- make_output_path(input_paths[[1]], ext, n_panels)
  ggsave(out_path, plot = combined,
         width = fig_w, height = fig_h, dpi = 300, bg = "white")
  message(sprintf("OK: %s", out_path))
}

tryCatch(
  main(input_paths, output_format),
  error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1)
  }
)

# =============================================================================
# JoMB STYLE SPECIFICATION (v9, applied to all panels)
# =============================================================================
# COLORS
#   #08306B navy    : mean fitness (density-scaled); beneficial mutations
#   #B1182C crimson : intrinsic fitness (genomic);   deleterious mutations
#
# TYPOGRAPHY
#   Family       : Arial
#   Axis titles  : bold, 15 pt
#   Axis ticks   : regular, 14 pt
#   Legend text  : regular, 14 pt
#   Panel tags   : bold, 19 pt (base_size + 4), top-left of each panel
#
# AXES & PANEL
#   Axis line     : black, 0.6 pt; ticks 4 pt long
#   No gridlines, no panel border
#   Y-axis: trajectory uses dynamic max; histogram uses log10 with shared max
#   X-axis: trajectory uses global max generation across panels
#
# LEGEND
#   Single shared legend at TOP of figure (via patchwork guides = "collect")
#   No border; centered
#
# LAYOUTS
#   2 files -> 1 x 2 side-by-side, 14 x 5.5 in
#   4 files -> 2 x 2 grid, 14 x 10 in
#
# OUTPUT
#   <basename_of_first_file>-<N>panel-<timestamp>.<ext>
#   300 dpi, white background
# =============================================================================
