#!/usr/bin/env Rscript
# slim_plot_generator.R
# -----------------------------------------------------------------------------
# Generates publication-grade plots from SLiM simulation outputs,
# styled to match the Rupe (JoMB) manuscript figures.
#
# Two file types are recognized by their first line:
#   1. Trajectory log    (CSV header begins with "cycle,Population_size")
#         -> Mean vs Intrinsic Fitness trajectory plot
#   2. Full output dump  (first line begins with "// Initial random seed:")
#         -> Allele frequency histogram (deleterious vs beneficial)
#
# Usage:
#   Rscript slim_plot_generator.R <file> [--format png|pdf|jpeg|tiff|svg] [--title "Plot title"]
#
# <path> must be a single .txt file. Output images are written next to the
# input as "<basename>-plot-<timestamp>.<format>".
#
# JoMB style summary (see bottom of file for full spec):
#   - Navy   #08306B : density-scaled mean fitness; beneficial mutations
#   - Crimson #B1182C : intrinsic genomic fitness;  deleterious mutations
#   - Times New Roman serif; bold axis titles; no gridlines
#   - Trajectory line width 0.4; legend top-LEFT (avoids overlap with rising
#     intrinsic fitness curves in the upper-right region)
#   - Histogram bin width 0.05; Fixed bin (>= 0.999) separated from the
#     0.0-1.0 range by a vertical dotted line at x = 1.0
# -----------------------------------------------------------------------------
# In R, `<-` is the conventional assignment operator. `=` works too, but `<-`
# is the idiom for binding values to names; `=` is reserved by convention for
# named function arguments. We follow that convention throughout.
# ---------- dependency loading ------------------------------------------------
# `suppressPackageStartupMessages` hides the chatter ggplot2/optparse print on
# library load - keeps the CLI output clean. `requireNamespace` returns FALSE
# instead of erroring if the package isn't installed, so we can give a friendly
# message before bailing.
require_pkg <- function(pkg) {
  # `requireNamespace(..., quietly = TRUE)` is the non-throwing "is it there?"
  # check. The `!` negates it: we enter the block when the package is missing.
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # `stop()` raises an error and (at the top level of an Rscript) terminates
    # with non-zero exit status. `call. = FALSE` suppresses the "Error in ..."
    # call-site prefix so the message reads cleanly.
    stop(sprintf("Missing R package '%s'. Install with: install.packages(\"%s\")",
                 pkg, pkg), call. = FALSE)
  }
}
require_pkg("optparse")   # CLI argument parser
require_pkg("ggplot2")    # plotting
require_pkg("scales")     # log-tick label formatter for histogram y-axis
# `suppressPackageStartupMessages({ ... })` wraps a block; libraries load
# silently. `library()` (vs `require()`) errors loudly if the package can't be
# attached - at this point we already verified it's installed, so any failure
# here is a real problem worth surfacing.
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(scales)
})
# ---------- JoMB style constants ----------------------------------------------
# These constants encode the visual style used throughout the Rupe (JoMB)
# manuscript figures. Centralizing them here makes it easy to tweak the look
# in one place rather than scattering hex codes and font sizes through the
# plotting code below.
JOMB_NAVY      <- "#08306B"   # mean fitness (density-scaled); beneficial muts
JOMB_CRIMSON   <- "#B1182C"   # intrinsic fitness (genomic);   deleterious muts
# Use Arial to match the table header font in the manuscript. On macOS,
# Arial is installed at /System/Library/Fonts/Supplemental/Arial.ttf and
# is available to R by name. If Arial isn't found, R falls back to the
# system sans-serif (Helvetica on macOS), which is visually near-identical.
JOMB_FONT      <- "Arial"
JOMB_BASE_SIZE <- 15          # base point size for axis titles (bold)

# Reusable theme that gives every plot the JoMB look in one call.
# Built on `theme_classic` (no gridlines, left and bottom axis lines only)
# with axis titles bolded, legend placed top-LEFT inside the panel with no
# border, and all text in Times New Roman serif. Top-left placement avoids
# overlap with intrinsic fitness curves that rise toward the upper-right
# (e.g., the H&C "Darwinian demon" trajectory).
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

      # Legend: top-LEFT inside panel, no border. Top-left avoids overlap
      # with rising intrinsic fitness curves in the upper-right region.
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
                                      margin = margin(b = 8))
    )
}
# ---------- CLI parsing -------------------------------------------------------
# `OptionParser` builds a parser much like argparse in Python or commander in
# JS. `make_option` declares one flag at a time. `type = "character"` means
# the value is captured as a string; `default` is used when the flag is absent.
option_list <- list(
  make_option(c("-f", "--format"),
              type = "character",
              default = "png",
              help = "Output image format: png (default), pdf, jpeg, tiff, svg",
              metavar = "FORMAT"),
  # --title is optional; when omitted (NULL) no title is rendered on the plot,
  # preserving the prior look. When provided, it's applied via `labs(title=...)`.
  make_option(c("-t", "--title"),
              type = "character",
              default = NULL,
              help = "Optional plot title (rendered above the chart). Default: no title.",
              metavar = "TITLE")
)
# `positional_arguments = 1` says: in addition to flags, expect exactly one
# positional argument (the input file path). `parse_args` returns a list
# with `$options` (named flags) and `$args` (the positional vector).
parser <- OptionParser(
  usage = "Usage: %prog <file> [--format png|pdf|jpeg|tiff|svg] [--title \"Plot title\"]",
  option_list = option_list
)
parsed <- tryCatch(
  parse_args(parser, positional_arguments = 1),
  # `error = function(e) ...` is R's tryCatch handler syntax. `e` is the
  # condition object; `conditionMessage(e)` extracts the human-readable text.
  error = function(e) {
    message("Argument error: ", conditionMessage(e))
    print_help(parser)
    quit(status = 2)   # `quit(status = N)` is R's process-exit-with-code.
  }
)
# Pull the positional argument and the named flags out of the parsed result.
input_path <- parsed$args[[1]]                # `[[1]]` extracts a single element by position
output_format <- tolower(parsed$options$format)  # normalize case for downstream comparisons
plot_title <- parsed$options$title            # NULL when --title was not supplied
# Whitelist the formats we know `ggsave()` can auto-detect from the extension.
# `%in%` is R's "value in vector" operator (like JS `arr.includes(x)`).
allowed_formats <- c("png", "pdf", "jpeg", "jpg", "tiff", "svg")
if (!(output_format %in% allowed_formats)) {
  stop(sprintf("Unsupported --format '%s'. Allowed: %s",
               output_format, paste(allowed_formats, collapse = ", ")),
       call. = FALSE)
}
# ---------- helpers -----------------------------------------------------------
# Detect file type by reading just the first line. Returns one of
# "trajectory", "fullout", or "unknown". Cheap (single line read) and the only
# heuristic we need given the two well-defined formats.
detect_type <- function(path) {
  # `readLines(path, n = 1)` reads up to one line. `warn = FALSE` silences the
  # "incomplete final line" warning R emits for files without a trailing \n.
  first_line <- tryCatch(
    readLines(path, n = 1, warn = FALSE),
    error = function(e) character(0)   # on read failure, return empty vector
  )
  # `length(first_line) == 0` means the file was empty or unreadable.
  if (length(first_line) == 0) return("unknown")
  # `startsWith(x, prefix)` is vectorized prefix-match (like JS String#startsWith).
  if (startsWith(first_line, "cycle,Population_size")) return("trajectory")
  if (startsWith(first_line, "// Initial random seed:")) return("fullout")
  return("unknown")
}
# Build the output path: same directory, basename + "-plot" + local timestamp
# + chosen extension. The timestamp avoids overwriting prior runs.
# `tools::file_path_sans_ext` strips the extension; `dirname` returns the parent
# directory. `file.path` joins path components with the OS separator.
# `format(Sys.time(), ...)` formats the current local-time POSIXct with strftime
# codes: %Y=year, %m=month, %d=day, %H=24h hour, %M=minute, %S=second.
make_output_path <- function(input_path, ext) {
  dir <- dirname(input_path)
  stem <- tools::file_path_sans_ext(basename(input_path))
  ts <- format(Sys.time(), "%Y_%m_%d-%H.%M.%S")   # local time, requested format
  file.path(dir, sprintf("%s-plot-%s.%s", stem, ts, ext))
}
# Resolve column names case-sensitively against the data frame, raising a
# clear error if any are missing. This is the centerpiece of the "no hardcoded
# indexes" rule for the trajectory CSV - columns are looked up by name only.
require_columns <- function(df, needed) {
  # `setdiff(A, B)` returns elements of A not in B - i.e. the missing ones.
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Trajectory file is missing required column(s): %s.\nColumns found: %s",
      paste(missing, collapse = ", "),
      paste(colnames(df), collapse = ", ")
    ), call. = FALSE)
  }
}
# ---------- Plot 1: trajectory ------------------------------------------------
plot_trajectory <- function(input_path, output_path, plot_title = NULL) {
  # `read.csv` with `check.names = FALSE` preserves header names verbatim
  # (otherwise R would mangle e.g. spaces or punctuation into dots).
  df <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  # Validate columns by NAME, never by position. SLiM users may add or
  # reorder columns; we only require these three exist.
  require_columns(df, c("cycle", "Mean_fitness", "Intrinsic_fitness"))
  # Build a "long" data frame: one row per (generation, series) pair. ggplot2
  # plots multi-series cleanly when the series identity is itself a column
  # (mapped via `aes(color = series)`), rather than two parallel columns.
  # We construct it manually with `rbind` of two `data.frame`s - no extra
  # package (tidyr/reshape2) needed.
  gen_kx <- df[["cycle"]] / 1000   # convert cycles to thousands of generations
  long_df <- rbind(
    data.frame(
      gen_kx = gen_kx,
      fitness = df[["Mean_fitness"]],
      series = "Mean fitness (density-scaled)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      gen_kx = gen_kx,
      fitness = df[["Intrinsic_fitness"]],
      series = "Intrinsic fitness (genomic)",
      stringsAsFactors = FALSE
    )
  )
  # `factor(..., levels = ...)` fixes the legend ordering (mean first, then
  # intrinsic) so colors map deterministically and the legend reads top-down
  # in our chosen order rather than alphabetically.
  long_df$series <- factor(
    long_df$series,
    levels = c("Mean fitness (density-scaled)", "Intrinsic fitness (genomic)")
  )
  # Manual color mapping: JoMB navy for mean, JoMB crimson for intrinsic.
  # Named vector keys must match the factor levels exactly.
  color_map <- c(
    "Mean fitness (density-scaled)" = JOMB_NAVY,
    "Intrinsic fitness (genomic)"   = JOMB_CRIMSON
  )
  # Build the ggplot object. In ggplot2, `+` chains layers/scales/themes onto
  # a base plot - analogous to method chaining in JS but additive.
  # Line width 0.4 matches the manuscript figures (was 0.3 in friend's original).
  p <- ggplot(long_df, aes(x = gen_kx, y = fitness, color = series)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_x_continuous(
      name = expression(bold(paste("Generation (", "\u00D7", "1,000)"))),
      limits = c(0, max(gen_kx, na.rm = TRUE)),   # `na.rm = TRUE` ignores NAs in max()
      expand = expansion(mult = c(0.005, 0.02))
    ) +
    scale_y_continuous(
      name = expression(bold("Fitness")),
      limits = c(0, 1.1),
      expand = expansion(mult = c(0, 0))
    ) +
    theme_jomb()
  # Apply optional title. `nzchar` rejects empty strings so `--title ""` doesn't
  # produce an empty title bar. `labs(title = ...)` is the standard ggplot2 way
  # to set/update the title without touching other labels.
  if (!is.null(plot_title) && nzchar(plot_title)) {
    p <- p + labs(title = plot_title)
  }
  # `ggsave` infers the device from the file extension. dpi only matters for
  # raster formats (png/jpeg/tiff); vector formats (pdf/svg) ignore it.
  # 8 x 6 inches matches the JoMB single-panel manuscript figure size.
  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
}
# ---------- Plot 2: allele frequency histogram --------------------------------
plot_allele_freq <- function(input_path, output_path, plot_title = NULL) {
  # Read the entire file into a character vector, one element per line.
  # SLiM full-output files are typically a few thousand to tens of thousands of
  # lines - tractable in memory.
  lines <- readLines(input_path, warn = FALSE)
  # ---- single-pass section walker ----
  # SLiM full-output sections appear in fixed order: Populations, Mutations,
  # Individuals, Genomes. We track which section we're currently inside via a
  # state variable. The section field positions below are POSITIONAL because
  # SLiM's spec is positional - see Reading Slim Full Output.md S3 for the
  # column meanings (mut_id at col 1, mut_type at col 3).
  section <- "header"
  # `mut_type_map` is a named character vector: names are mut_ids (as strings),
  # values are "m2"/"m3"/etc. Pre-allocating the right size would be faster,
  # but for typical full-output files (a few thousand mutations) the simple
  # accumulator is fine. We collect into lists and convert at the end to avoid
  # quadratic-cost growth-by-append.
  mut_ids_acc   <- character(0)   # accumulator for mut_id strings
  mut_types_acc <- character(0)   # parallel accumulator for mut_type strings
  # Genome parsing accumulates ALL referenced mut_ids across all genome lines
  # into one big vector, then `table()` counts in a single pass at the end.
  # This avoids the O(N^2) cost of incrementing a named-vector counter inside
  # a loop (R's named-vector index updates are not in-place).
  genome_count <- 0L
  genome_ids_chunks <- list()   # list of integer vectors, one per genome line
  # Pre-compile the "starts with 'p1:'" check; `startsWith` is already fast,
  # but we hoist the constant out of the loop for readability.
  for (raw in lines) {
    # Section headers in SLiM full output are bare lines like "Mutations:".
    # Detect by exact match; cheaper than regex.
    if (raw == "Mutations:")        { section <- "mutations";   next }
    if (raw == "Individuals:")      { section <- "individuals"; next }
    if (raw == "Genomes:")          { section <- "genomes";     next }
    if (startsWith(raw, "Populations:")) { section <- "populations"; next }
    # Skip empty lines and the header preamble (// comments, #OUT, Version).
    if (section == "header" || section == "populations" || section == "individuals") next
    if (!nzchar(raw)) next   # `nzchar` = "non-zero character count"
    if (section == "mutations") {
      # `strsplit(x, "\\s+")` splits on runs of whitespace and returns a list;
      # `[[1]]` unwraps the single result (since we pass one string in).
      fields <- strsplit(raw, "\\s+")[[1]]
      # SLiM Mutations spec (positional, per format doc):
      #   fields[1] = file-local mut_id (referenced in Genomes section)
      #   fields[3] = mutation type (e.g. "m2", "m3")
      # Other fields (permanent ID, position, s, h, origin subpop, origin tick,
      # prevalence) exist but are not needed for this plot.
      if (length(fields) >= 3) {
        mut_ids_acc[length(mut_ids_acc) + 1L]     <- fields[1]
        mut_types_acc[length(mut_types_acc) + 1L] <- fields[3]
      }
      next
    }
    if (section == "genomes") {
      # Genome lines look like: "p1:0 A 0 1 2 3 ... 509"
      # Skip any malformed lines that don't start with "p1:".
      if (!startsWith(raw, "p1:")) next
      genome_count <- genome_count + 1L
      fields <- strsplit(raw, "\\s+")[[1]]
      # fields[1] = "p1:N" (genome ID); fields[2] = chromosome type ("A"/"X"/"Y").
      # Remaining fields (if any) are mut_ids referenced by this genome.
      if (length(fields) > 2) {
        # `as.integer` vectorizes over the slice; produces NA on parse failure
        # which we strip with `na.omit`-equivalent indexing.
        ids <- suppressWarnings(as.integer(fields[-c(1, 2)]))
        ids <- ids[!is.na(ids)]
        # Append the chunk to our list. Growing a list by index is O(1) amortized
        # in R (unlike named-vector growth); we'll concatenate once at the end.
        genome_ids_chunks[[length(genome_ids_chunks) + 1L]] <- ids
      }
      next
    }
  }
  # ---- post-loop consolidation ----
  # Sanity check: refuse to plot a degenerate file.
  if (genome_count == 0L) {
    stop(sprintf("No genome lines found in '%s' - is this a full-output file?",
                 input_path), call. = FALSE)
  }
  if (length(mut_ids_acc) == 0L) {
    stop(sprintf("No mutations parsed from '%s'.", input_path), call. = FALSE)
  }
  # Build the mut_id -> mut_type lookup. Names of the vector ARE the mut_ids;
  # `mut_type_map["167"]` returns "m2" etc.
  mut_type_map <- mut_types_acc
  names(mut_type_map) <- mut_ids_acc
  # Concatenate all genome chunks into one big integer vector, then count.
  # `unlist` flattens a list of vectors; `table` produces a contingency count
  # (named integer vector: name = mut_id as string, value = how many genomes
  # carry it). This is the efficient path mentioned in the design.
  all_ids <- unlist(genome_ids_chunks, use.names = FALSE)
  counts <- table(all_ids)   # named integer vector keyed by mut_id (as string)
  # Frequency = count / total_genomes. Element-wise division is vectorized.
  freqs <- as.numeric(counts) / genome_count
  freq_ids <- names(counts)   # mut_ids as strings, parallel to `freqs`
  # Look up each mut_id's type. `mut_type_map[freq_ids]` does a vectorized
  # name-based lookup; mut_ids referenced by genomes but missing from the
  # Mutations section (shouldn't happen in a well-formed file) yield NA.
  types_for_freq <- mut_type_map[freq_ids]
  # Split into deleterious (m2) and beneficial (m3) frequency vectors.
  # Other types (m1 neutral, etc.) are ignored per the pseudocode spec.
  del_freqs <- freqs[!is.na(types_for_freq) & types_for_freq == "m2"]
  ben_freqs <- freqs[!is.na(types_for_freq) & types_for_freq == "m3"]
  # ---- JoMB binning spec ----
  # Fixed mutations: frequency >= 0.999 (essentially 1.0).
  # Non-fixed mutations: binned in 0.05-wide bins from 0.00 to 1.00.
  # The Fixed bin is plotted separately at a position to the right of 1.0,
  # with a vertical dotted line at x = 1.0 acting as a visual separator.
  FIXED_THRESHOLD <- 0.999
  BIN_WIDTH       <- 0.05
  # Fixation counts go into the legend labels.
  del_fixed <- sum(del_freqs >= FIXED_THRESHOLD)
  ben_fixed <- sum(ben_freqs >= FIXED_THRESHOLD)
  # Separate fixed from non-fixed for binning.
  del_nonfixed <- del_freqs[del_freqs < FIXED_THRESHOLD]
  ben_nonfixed <- ben_freqs[ben_freqs < FIXED_THRESHOLD]
  # Build the bin edges. `seq(0, 1, by = 0.05)` gives 21 edges (0.00, 0.05,
  # 0.10, ..., 1.00), which is 20 bins. Edges are right-inclusive on the last
  # bin via `include.lowest = TRUE` below (so a frequency of exactly 0 lands
  # in the first bin).
  edges <- seq(0, 1, by = BIN_WIDTH)
  # Bin centers (used as continuous x-coordinates so we can place the Fixed
  # bin to the right of x = 1 with a real numeric gap).
  bin_centers <- edges[-length(edges)] + BIN_WIDTH / 2
  # `cut` returns a factor of bin assignments. We pass `labels = FALSE` to
  # get integer bin indices instead - more convenient for counting via tabulate.
  bin_idx_del <- cut(del_nonfixed, breaks = edges,
                     include.lowest = TRUE, right = FALSE, labels = FALSE)
  bin_idx_ben <- cut(ben_nonfixed, breaks = edges,
                     include.lowest = TRUE, right = FALSE, labels = FALSE)
  # `tabulate(x, nbins = N)` returns counts for each integer from 1..N,
  # including zero counts for empty bins. This gives parallel vectors we can
  # zip into the long-format plot data frame.
  n_bins <- length(bin_centers)
  del_counts <- tabulate(bin_idx_del, nbins = n_bins)
  ben_counts <- tabulate(bin_idx_ben, nbins = n_bins)
  # ---- assemble plotting data frame ----
  # Long form: one row per (bin_center, category) cell. Replace zero counts
  # with NA so the log-y scale doesn't emit "removed N rows containing missing
  # values" warnings AND so empty bars simply don't render (cleaner than a
  # -Inf bar).
  zero_to_na <- function(v) ifelse(v == 0L, NA_integer_, v)
  # Legend labels include the fixed counts via `sprintf` (R's printf).
  del_label <- sprintf("Deleterious (fixed: %d)", del_fixed)
  ben_label <- sprintf("Beneficial (fixed: %d)", ben_fixed)
  # Position the Fixed bin two bin-widths to the right of x = 1.0 so there's
  # a clear visual gap. The vertical dotted separator (added below as a
  # geom_vline) sits between x = 1.0 and the Fixed bin.
  FIXED_X <- 1.0 + 2 * BIN_WIDTH
  # Build the data frame: numeric bin centers for non-fixed bars, plus one
  # row per category for the Fixed bar at FIXED_X.
  plot_df <- rbind(
    data.frame(
      bin_x = bin_centers,
      count = zero_to_na(del_counts),
      category = del_label,
      stringsAsFactors = FALSE
    ),
    data.frame(
      bin_x = bin_centers,
      count = zero_to_na(ben_counts),
      category = ben_label,
      stringsAsFactors = FALSE
    ),
    data.frame(
      bin_x = FIXED_X,
      count = ifelse(del_fixed == 0L, NA_integer_, del_fixed),
      category = del_label,
      stringsAsFactors = FALSE
    ),
    data.frame(
      bin_x = FIXED_X,
      count = ifelse(ben_fixed == 0L, NA_integer_, ben_fixed),
      category = ben_label,
      stringsAsFactors = FALSE
    )
  )
  # Lock category order so deleterious appears first (crimson, on the left of
  # each dodge group), beneficial second (navy on the right).
  plot_df$category <- factor(plot_df$category,
                             levels = c(del_label, ben_label))
  # Color map keyed by the dynamic legend labels (which include fix counts).
  # JoMB crimson for deleterious, JoMB navy for beneficial.
  fill_map <- setNames(c(JOMB_CRIMSON, JOMB_NAVY), c(del_label, ben_label))
  # Bar width: slightly narrower than the bin width so adjacent bins don't
  # touch. Dodge width matches bin width so the deleterious/beneficial pair
  # within a bin sits centered on the bin center.
  # Bars sized to be thick and visible: each bar = 80% of bin width, with
  # dodge equal to bar width so the two bars in a pair sit immediately
  # adjacent. Because bar_width and dodge_width are equal, the deleterious
  # and beneficial bars within each bin sit flush against each other.
  bar_width   <- BIN_WIDTH * 0.8
  dodge_width <- BIN_WIDTH * 0.8
  # X-axis ticks: every 0.1 from 0.0 to 1.0, plus a "Fixed" tick at FIXED_X.
  tick_positions <- c(seq(0, 1, by = 0.1), FIXED_X)
  tick_labels    <- c(sprintf("%.1f", seq(0, 1, by = 0.1)), "Fixed")
  # Build the plot. `geom_col` draws bars at the data's y values (vs `geom_bar`
  # which counts rows). `position_dodge` places paired bars side-by-side.
  p <- ggplot(plot_df, aes(x = bin_x, y = count, fill = category)) +
    geom_col(position = position_dodge(width = dodge_width),
             width = bar_width, na.rm = TRUE, color = NA) +
    # Vertical dotted separator between the numeric bins and the Fixed bin.
    geom_vline(xintercept = 1.0 + BIN_WIDTH / 2,
               linetype = "dotted", color = "black", linewidth = 0.4) +
    scale_fill_manual(values = fill_map, name = NULL) +
    scale_y_log10(
      name = expression(bold("Number of mutations (log scale)")),
      labels = trans_format("log10", math_format(10^.x)),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(
      name = expression(bold("Allele frequency")),
      breaks = tick_positions,
      labels = tick_labels,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    theme_jomb() +
    # Histogram-only legend tweak: nudge the legend rightward (away from the
    # tall low-frequency bars) and use the same default ggplot legend keys
    # as the trajectory plot so font, size, and spacing match.
    theme(legend.position = c(0.12, 0.97))
  # Apply optional title (same idiom as in plot_trajectory).
  if (!is.null(plot_title) && nzchar(plot_title)) {
    p <- p + labs(title = plot_title)
  }
  # Histograms are slightly wider than trajectories to give room for the
  # Fixed bin and the 0.1-step x-axis labels.
  ggsave(output_path, plot = p, width = 9, height = 6, dpi = 300, bg = "white")
}
# ---------- dispatch ----------------------------------------------------------
# Top-level entry point. Single-file only (folder mode was removed by request).
# Steps: existence check -> reject directories -> require .txt extension ->
# detect format from first line -> dispatch to the matching plot function ->
# save. Any failure raises an error that the outer tryCatch turns into a clean
# non-zero exit.
main <- function(input_path, ext, plot_title) {
  # `file.exists` is true for both files and directories; combined with the
  # `dir.exists` check below it lets us reject the directory case explicitly.
  if (!file.exists(input_path)) {
    stop(sprintf("Path does not exist: %s", input_path), call. = FALSE)
  }
  if (dir.exists(input_path)) {
    stop(sprintf(
      "Folder input is not supported - pass a single .txt file. Got: %s",
      input_path
    ), call. = FALSE)
  }
  # Per spec: if the extension isn't .txt, error out without even reading the
  # file. `tools::file_ext` returns the extension sans dot; lower-cased compare.
  if (tolower(tools::file_ext(input_path)) != "txt") {
    stop(sprintf("Input file must have .txt extension: %s", input_path),
         call. = FALSE)
  }
  kind <- detect_type(input_path)
  if (kind == "unknown") {
    stop(sprintf(
      "Unrecognized SLiM file format (first line matched neither trajectory nor full-output): %s",
      input_path
    ), call. = FALSE)
  }
  out <- make_output_path(input_path, ext)
  if (kind == "trajectory") {
    plot_trajectory(input_path, out, plot_title = plot_title)
  } else {
    plot_allele_freq(input_path, out, plot_title = plot_title)
  }
  message(sprintf("OK: %s", out))
}
# Run main with top-level error handling. Any uncaught error from `main`
# terminates the script with a clean message and exit code 1.
tryCatch(
  main(input_path, output_format, plot_title),
  error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1)
  }
)
# =============================================================================
# JoMB STYLE SPECIFICATION (applied throughout this script)
# =============================================================================
#
# COLORS
#   #08306B navy    : mean fitness (density-scaled); beneficial mutations
#   #B1182C crimson : intrinsic fitness (genomic);   deleterious mutations
#
# TYPOGRAPHY
#   Family       : Arial (matches manuscript table header font)
#   Axis titles  : Arial bold, 15 pt
#   Axis ticks   : Arial regular, 14 pt, black
#   Legend text  : Arial regular, 14 pt, no title
#
# AXES & PANEL
#   Axis line    : black, 0.6 pt; ticks 4 pt long
#   No gridlines, no panel border
#   Y axis starts at 0 with no padding below
#
# LEGEND
#   Position : top-LEFT inside panel (0.03, 0.97 with left/top justification)
#   No border (matches manuscript figures); white background
#   Top-LEFT (not top-right) so legend doesn't overlap with rising intrinsic
#   fitness curves in the upper-right region (e.g., the Darwinian demon case).
#
# TRAJECTORY PLOTS
#   Line width: 0.7 pt
#   X-axis label: "Generation (x1,000)"
#   Y-axis label: "Fitness"
#   Y-axis range: [0, 1.1]
#   Output size : 8 in x 6 in at 300 dpi
#
# HISTOGRAMS
#   Bin width: 0.05 (20 bins across [0, 1])
#   Fixed bin: frequencies >= 0.999, placed at synthetic x = 1.0 + 2 * bin_width
#   Separator: vertical dotted black line at x = 1.0 + bin_width/2
#   Dodged bars, no outline; bar width = 0.4 x bin_width = 0.02
#   Y-axis: log10 with 10^N tick labels
#   X-axis ticks: 0.0, 0.1, ..., 1.0, Fixed
#   Output size : 9 in x 6 in at 300 dpi
#
# OUTPUT
#   Background: white (set via bg = "white" in ggsave to avoid transparency)
#
# =============================================================================
