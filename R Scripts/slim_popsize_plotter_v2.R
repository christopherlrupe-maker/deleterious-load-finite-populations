#!/usr/bin/env Rscript
# =============================================================================
# slim_popsize_plotter_v2.R
# =============================================================================
# Standalone command-line tool that reads a SLiM trajectory log (.txt) and
# produces a population-size-over-time plot matching the manuscript style
# (v14 styling: Arial bold axis titles, navy line, top-left vertical legend,
# no border).
#
# Expected input columns (case-sensitive):
#   cycle             -- generation number
#   Population_size   -- population size at that generation
#
# Usage:
#   Rscript slim_popsize_plotter_v2.R "<path-to-file.txt>" \
#           [--title "Plot title"] \
#           [--format png|pdf|jpeg|tiff|svg] \
#           [--k 1000]
#
# Output:
#   Saved next to the input file as:
#     <original-name>-popsize-YYYY_MM_DD-HH.MM.SS.<format>
# =============================================================================

# ---------- package loading ---------------------------------------------------
require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Missing R package '%s'. Install with: install.packages(\"%s\")",
                 pkg, pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
require_pkg("optparse")
require_pkg("ggplot2")
require_pkg("scales")

# ---------- JoMB style constants (matching v14) -------------------------------
JOMB_NAVY      <- "#08306B"
JOMB_CRIMSON   <- "#B1182C"
JOMB_GRAY      <- "gray55"    # for dashed reference lines
JOMB_FONT      <- "Arial"
JOMB_BASE_SIZE <- 15
JOMB_LINE_W    <- 0.35        # trajectory line width

# Reusable theme that gives every plot the JoMB look (same as v14).
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

      # Legend: top-left inside panel, no border (matches v14 trajectory plot)
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
      plot.title       = element_text(hjust = 0.5, face = "bold",
                                      size = base_size + 1,
                                      margin = margin(b = 8))
    )
}

# ---------- argument parsing --------------------------------------------------
option_list <- list(
  make_option("--title",  type = "character", default = NULL,
              help = "Optional plot title"),
  make_option("--format", type = "character", default = "png",
              help = "Output image format: png, pdf, jpeg, tiff, svg [default %default]"),
  make_option("--k",      type = "double",    default = 1000,
              help = "Carrying capacity for dashed reference line [default %default]")
)
parser <- OptionParser(usage = "Rscript slim_popsize_plotter_v2.R <input.txt> [options]",
                       option_list = option_list)
args <- parse_args(parser, positional_arguments = 1)
input_path <- args$args[1]
opts       <- args$options

# ---------- input validation --------------------------------------------------
if (!file.exists(input_path)) {
  stop(sprintf("Input file not found: %s", input_path), call. = FALSE)
}
if (dir.exists(input_path)) {
  stop("Folder input is not supported. Pass a single .txt file.", call. = FALSE)
}
if (!grepl("\\.txt$", input_path, ignore.case = TRUE)) {
  stop("Input file must have .txt extension.", call. = FALSE)
}

valid_formats <- c("png", "pdf", "jpeg", "tiff", "svg")
if (!(tolower(opts$format) %in% valid_formats)) {
  stop(sprintf("Unrecognized --format '%s'. Choose one of: %s",
               opts$format, paste(valid_formats, collapse = ", ")),
       call. = FALSE)
}

# ---------- read and validate data --------------------------------------------
df <- tryCatch(
  read.csv(input_path, stringsAsFactors = FALSE),
  error = function(e) stop(sprintf("Could not read file as CSV: %s", e$message),
                           call. = FALSE)
)

needed <- c("cycle", "Population_size")
missing_cols <- setdiff(needed, names(df))
if (length(missing_cols) > 0) {
  stop(sprintf("Trajectory file is missing required column(s): %s. Columns found: %s",
               paste(missing_cols, collapse = ", "),
               paste(names(df), collapse = ", ")),
       call. = FALSE)
}

# ---------- build plot --------------------------------------------------------
# Convert generation to thousands so the x-axis matches the other figures.
df$gen_k <- df$cycle / 1000

# Legend label for carrying capacity reference line. format() inserts commas.
k_label   <- sprintf("Carrying capacity (K = %s)", format(opts$k, big.mark = ","))

# Use aes(linetype = ...) on geom_hline so the dashed reference line appears
# in the legend alongside the solid population-size line. We map "Population
# size" to the color aesthetic and the K reference to the linetype aesthetic;
# both legends are collected into a single visual block.
p <- ggplot(df, aes(x = gen_k, y = Population_size)) +
  geom_hline(aes(yintercept = opts$k, linetype = k_label),
             color = JOMB_GRAY, linewidth = 0.4) +
  geom_line(aes(color = "Population size"), linewidth = JOMB_LINE_W) +
  scale_color_manual(name = NULL,
                     values = c("Population size" = JOMB_NAVY)) +
  scale_linetype_manual(name = NULL,
                        values = stats::setNames("dashed", k_label)) +
  scale_x_continuous(
    name = expression(bold(paste("Generation (", "\u00D7", "1,000)"))),
    expand = expansion(mult = c(0.005, 0.02))
  ) +
  scale_y_continuous(
    name = expression(bold("Population size")),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_jomb() +
  # Force consistent legend rendering: both keys in one block, navy line
  # at full width for the population trace, dashed gray line for K.
  guides(
    color    = guide_legend(order = 1,
                            override.aes = list(linewidth = JOMB_LINE_W * 2)),
    linetype = guide_legend(order = 2,
                            override.aes = list(color = JOMB_GRAY,
                                                linewidth = 0.4))
  )

# Apply optional title (same idiom as v14)
if (!is.null(opts$title) && nzchar(opts$title)) {
  p <- p + labs(title = opts$title)
}

# ---------- save --------------------------------------------------------------
stamp     <- format(Sys.time(), "%Y_%m_%d-%H.%M.%S")
base_name <- sub("\\.txt$", "", basename(input_path), ignore.case = TRUE)
out_dir   <- dirname(input_path)
out_name  <- sprintf("%s-popsize-%s.%s", base_name, stamp, tolower(opts$format))
out_path  <- file.path(out_dir, out_name)

# 8 x 6 in at 300 dpi (same as v14 trajectory output)
ggsave(out_path, plot = p, width = 8, height = 6, dpi = 300, bg = "white")

cat(sprintf("OK: %s\n", out_path))

# =============================================================================
# STYLE NOTES (matches v14)
# =============================================================================
#   Navy        #08306B  : main population-size line
#   Gray55      gray55   : dashed carrying-capacity reference line
#   Font        Arial    : axis titles bold 15 pt, ticks/legend 14 pt regular
#   Legend      top-left (c(0.03, 0.97)), no border
#   Line width  0.35 pt  : trajectory line
#   Output      8 x 6 in at 300 dpi, white background
# =============================================================================
