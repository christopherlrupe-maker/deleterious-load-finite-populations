# =============================================================================
# DFE of fixed deleterious mutations during the stable phase — JoMB Styling
# =============================================================================
# Four-panel histogram comparing the distribution of |s| for fixed deleterious
# mutations in:
#   a) H&C baseline (5 replicates pooled, mean s = -0.5)
#   b) B1 Run 1 stable phase (mean s = -0.1)
#   c) B1 Run 2 stable phase (mean s = -0.05)
#   d) B1 Run 3 stable phase (mean s = -0.01)
#
# Each panel shades the effectively neutral zone (ENZ, |s| < 0.0005, pale green)
# and the selection breakdown zone (SBZ, 0.0005 <= |s| <= 0.0025, pale orange),
# matching the styling of DFE_comparison_ENZ_SBZ.R (v14 JoMB family).
#
# Style notes:
#   - Arial font, axis titles bold 15 pt, ticks 14 pt
#   - Axis lines black 0.6 pt, 4 pt ticks, no gridlines
#   - 10 x 8 in at 300 dpi, white background (2x2 grid)
#   - Per-panel histogram color matches the run color from the comparison plot
#     (red / blue / green / purple)
#   - Panel tags (a-d) bold Arial, outside upper-left of each panel
#   - ENZ + SBZ zone legend appears in panel (a) only
#
# Methodology
# -----------
# H&C baseline: pooled tab-separated mutation table (GE_run_full.txt).
#   Keep rows with mut_type == "m2" (deleterious) and frequency == 1 (fixed).
#   Pool across all 5 replicates. The entire 50k-gen run is stable phase, so
#   no gen_origin filter. |s| = |effect|.
#
# B1 runs: raw SLiM output (single near-extinction dump per run). The
#   "Mutations:" section records segregating mutations as whitespace-separated:
#     tempID mutID mut_type pos s_coef h subpop gen_origin numCopies
#   At dump time the population is N = 4 (2N = 8 chromosomes), so a mutation
#   with numCopies == 8 is fixed in the surviving population. Keep m2 fixations
#   whose gen_origin precedes the stable-phase cutoff (the last generation
#   before intrinsic fitness dropped below 0.1, derived empirically from the
#   per-cycle intrinsic-fitness traces):
#       B1R1 (mean s = -0.1):   gen_origin <= 203500   (~92% of the run)
#       B1R2 (mean s = -0.05):  gen_origin <=  66300   (~83% of the run)
#       B1R3 (mean s = -0.01):  gen_origin <=  30900   (~79% of the run)
#
# Requires: ggplot2, patchwork
#   install.packages(c("ggplot2", "patchwork"))
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})


# ---------- JoMB style constants (matching v14) ------------------------------
JOMB_FONT      <- "Arial"
JOMB_BASE_SIZE <- 15

# Fall back to a generic sans-serif if Arial isn't installed locally (e.g.,
# on a sandboxed CI runner). On the author's machine where Arial IS
# installed, this leaves JOMB_FONT unchanged.
.has_arial <- tryCatch(
  {
    sf <- systemfonts::match_fonts("Arial")
    isTRUE(grepl("Arial", sf$family, ignore.case = TRUE))
  },
  error = function(e) FALSE
)
if (!.has_arial) {
  message("Arial not detected; falling back to 'sans' for this environment.")
  JOMB_FONT <- "sans"
}

theme_jomb <- function(base_size = JOMB_BASE_SIZE,
                       base_family = JOMB_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Axis titles slightly smaller than the base size so they sit lighter
      # alongside the panel tags.
      axis.title.x = element_text(face = "bold", size = base_size - 2,
                                  margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = base_size - 2,
                                  margin = margin(r = 8)),
      axis.text.x  = element_text(size = base_size - 2, color = "black"),
      axis.text.y  = element_text(size = base_size - 2, color = "black"),
      axis.line    = element_line(color = "black", linewidth = 0.6),
      axis.ticks   = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(4, "pt"),

      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),

      legend.background = element_blank(),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size - 2,
                                       family = base_family),

      plot.background = element_rect(fill = "white", color = NA),
      # Top margin needs to be generous because cowplot places each panel's
      # label (a/b/c/d) in this space, just outside the upper-left corner.
      # Right margin needs extra room so the rightmost x-tick label ("0.005")
      # isn't clipped by the panel boundary.
      plot.margin     = margin(t = 30, r = 22, b = 8, l = 14)
    )
}


# ---------- Configuration ----------------------------------------------------
# This script supports two ways of telling it where the four input files
# live. The four files it needs are:
#
#   --hc    GE_run_full.txt                                (H&C baseline)
#   --b1r1  GE_run_mean_fit_-0_1_300k_Full_Output.txt      (mean s = -0.1)
#   --b1r2  GE_run_mean_fit_-0_05_100k_Full_Output.txt     (mean s = -0.05)
#   --b1r3  GE_run_mean_-0_01_50k_Full_Output.txt          (mean s = -0.01)
#
# (A) Single-directory mode — all four files in the same folder:
#       Rscript DFE_fixed_deleterious_JoMB.R /path/to/data
#       Rscript DFE_fixed_deleterious_JoMB.R /path/to/data /path/to/output
#     Output defaults to the input directory.
#
# (B) Per-file mode — files in different folders (use named flags):
#       Rscript DFE_fixed_deleterious_JoMB.R \
#         --hc   /path/.../HC_panel/GE_run_full.txt \
#         --b1r1 /path/.../B1R1_panel/GE_run_mean_fit_-0_1_300k_Full_Output.txt \
#         --b1r2 /path/.../B1R2_panel/GE_run_mean_fit_-0_05_100k_Full_Output.txt \
#         --b1r3 /path/.../B1R3_panel/GE_run_mean_-0_01_50k_Full_Output.txt \
#         --out  /path/to/output/dir
#     If --out is omitted, output is written to the directory of --hc.

.cli_args <- commandArgs(trailingOnly = TRUE)

# Helper: pull "--flag VALUE" out of an argument vector
.flag_value <- function(args, flag) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) NA_character_ else args[[i + 1L]]
}

.flags <- c("--hc", "--b1r1", "--b1r2", "--b1r3", "--out")
.has_flags <- any(.flags %in% .cli_args)

FILE_PATHS <- list()
if (.has_flags) {
  # Per-file mode (B)
  FILE_PATHS$hc   <- .flag_value(.cli_args, "--hc")
  FILE_PATHS$b1r1 <- .flag_value(.cli_args, "--b1r1")
  FILE_PATHS$b1r2 <- .flag_value(.cli_args, "--b1r2")
  FILE_PATHS$b1r3 <- .flag_value(.cli_args, "--b1r3")
  OUTPUT_DIR <- .flag_value(.cli_args, "--out")
  if (is.na(OUTPUT_DIR) || !nzchar(OUTPUT_DIR)) {
    OUTPUT_DIR <- if (!is.na(FILE_PATHS$hc)) dirname(FILE_PATHS$hc) else getwd()
  }
} else {
  # Single-directory mode (A)
  INPUT_DIR <- if (length(.cli_args) >= 1L) {
    .cli_args[[1L]]
  } else {
    script_path <- tryCatch(
      {
        cmd_args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", cmd_args, value = TRUE)
        if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
      },
      error = function(e) NA_character_
    )
    if (!is.na(script_path) && nzchar(script_path)) {
      dirname(normalizePath(script_path, mustWork = FALSE))
    } else {
      getwd()
    }
  }
  OUTPUT_DIR <- if (length(.cli_args) >= 2L) .cli_args[[2L]] else INPUT_DIR
  FILE_PATHS$hc   <- file.path(INPUT_DIR, "GE_run_full.txt")
  FILE_PATHS$b1r1 <- file.path(INPUT_DIR, "GE_run_mean_fit_-0_1_300k_Full_Output.txt")
  FILE_PATHS$b1r2 <- file.path(INPUT_DIR, "GE_run_mean_fit_-0_05_100k_Full_Output.txt")
  FILE_PATHS$b1r3 <- file.path(INPUT_DIR, "GE_run_mean_-0_01_50k_Full_Output.txt")
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

OUTPUT_PNG <- file.path(OUTPUT_DIR, "DFE_fixed_deleterious.png")
OUTPUT_PDF <- file.path(OUTPUT_DIR, "DFE_fixed_deleterious.pdf")

message("Input files:")
for (nm in names(FILE_PATHS)) {
  message(sprintf("  %-5s %s", paste0(nm, ":"), FILE_PATHS[[nm]]))
}
message("Output PNG:  ", OUTPUT_PNG)
message("Output PDF:  ", OUTPUT_PDF)

ENZ_MAX <- 0.0005
SBZ_MAX <- 0.0025
X_MAX   <- 0.005
N_BINS  <- 50
BIN_EDGES <- seq(0, X_MAX, length.out = N_BINS + 1)

# Population is N = 4 individuals (2N = 8 chromosomes) at the dump moment
FIXED_COPIES <- 8

# Stable-phase cutoffs: last generation before intrinsic fitness < 0.1
CUTOFFS <- list(
  B1R1 = 203500,
  B1R2 =  66300,
  B1R3 =  30900
)

# Run colors — match the v14 comparison-plot palette so panels are
# immediately recognizable across figures in the manuscript.
RUN_COLORS <- list(
  HC   = "#C9302C",   # red
  B1R1 = "#1F4E79",   # blue
  B1R2 = "#1F6F2F",   # green
  B1R3 = "#5E2D8A"    # purple
)

# Zone shading (matches v14)
COL_ENZ_FILL <- "#C8E6B8"   # pale green
COL_SBZ_FILL <- "#FFCC99"   # pale orange
COL_ENZ_LINE <- "#2E7D32"   # dark green
COL_SBZ_LINE <- "#E07000"   # dark orange


# ---------- Helper: read H&C pooled mutation table ---------------------------
read_hc_fixed_s <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  fixed_del <- df[df$mut_type == "m2" & df$frequency == 1, ]
  abs(fixed_del$effect)
}


# ---------- Helper: extract fixed m2 in stable phase from a SLiM dump --------
read_b1_fixed_s <- function(path, gen_cutoff) {
  lines <- readLines(path)
  mut_start <- which(lines == "Mutations:") + 1L
  mut_end   <- which(startsWith(lines, "Individuals:")) - 1L
  stopifnot(length(mut_start) == 1, length(mut_end) == 1)

  fields <- strsplit(trimws(lines[mut_start:mut_end]), "\\s+")
  ok <- vapply(fields, length, integer(1)) >= 9
  fields <- fields[ok]

  mut_type   <- vapply(fields, `[`, character(1), 3L)
  s_coef     <- as.numeric(vapply(fields, `[`, character(1), 5L))
  gen_origin <- as.integer(vapply(fields, `[`, character(1), 8L))
  num_copies <- as.integer(vapply(fields, function(x) x[length(x)], character(1)))

  keep <- mut_type == "m2" &
          num_copies == FIXED_COPIES &
          gen_origin <= gen_cutoff
  abs(s_coef[keep])
}


# ---------- Load data --------------------------------------------------------
# Verify all four files exist; if any is missing, print a clear error before
# attempting to read.
.missing <- vapply(FILE_PATHS, function(p) is.na(p) || !file.exists(p), logical(1))
if (any(.missing)) {
  msg <- "Missing or unspecified input file(s):\n"
  for (nm in names(FILE_PATHS)[.missing]) {
    msg <- paste0(msg, sprintf("  %-5s %s\n", paste0(nm, ":"),
                               if (is.na(FILE_PATHS[[nm]])) "(not provided)"
                               else FILE_PATHS[[nm]]))
  }
  msg <- paste0(
    msg,
    "\nProvide the four files via one of:\n",
    "  Rscript DFE_fixed_deleterious_JoMB.R /path/to/folder/with/all/four\n",
    "  Rscript DFE_fixed_deleterious_JoMB.R \\\n",
    "    --hc   /path/to/GE_run_full.txt \\\n",
    "    --b1r1 /path/to/GE_run_mean_fit_-0_1_300k_Full_Output.txt \\\n",
    "    --b1r2 /path/to/GE_run_mean_fit_-0_05_100k_Full_Output.txt \\\n",
    "    --b1r3 /path/to/GE_run_mean_-0_01_50k_Full_Output.txt \\\n",
    "    --out  /path/to/output/dir"
  )
  stop(msg, call. = FALSE)
}

hc   <- read_hc_fixed_s(FILE_PATHS$hc)
b1r1 <- read_b1_fixed_s(FILE_PATHS$b1r1, CUTOFFS$B1R1)
b1r2 <- read_b1_fixed_s(FILE_PATHS$b1r2, CUTOFFS$B1R2)
b1r3 <- read_b1_fixed_s(FILE_PATHS$b1r3, CUTOFFS$B1R3)

# Sanity check vs. the target figure
for (nm in c("hc", "b1r1", "b1r2", "b1r3")) {
  v <- get(nm)
  pct <- 100 * sum(v <= SBZ_MAX) / length(v)
  message(sprintf("%-5s: n = %4d, ENZ+SBZ = %.1f%%", nm, length(v), pct))
}


# ---------- Per-panel metadata -----------------------------------------------
fmt_count <- function(n) format(n, big.mark = ",")

panels <- list(
  list(label = "a",
       data  = hc,
       color = RUN_COLORS$HC,
       y_max =  600, y_step = 100,
       box = c(
         "H&C baseline (5 replicates pooled)",
         "(mean s = \u20130.5)",
         paste0(fmt_count(length(hc)), " fixed deleterious"),
         sprintf("%.1f%% in ENZ + SBZ", 100 * sum(hc <= SBZ_MAX) / length(hc))
       ),
       show_zone_legend = TRUE),

  list(label = "b",
       data  = b1r1,
       color = RUN_COLORS$B1R1,
       y_max = 1200, y_step = 200,
       box = c(
         "B1 Run 1 stable phase",
         "(mean s = \u20130.1)",
         paste0(fmt_count(length(b1r1)), " fixed deleterious"),
         sprintf("%.1f%% in ENZ + SBZ", 100 * sum(b1r1 <= SBZ_MAX) / length(b1r1))
       ),
       show_zone_legend = FALSE),

  list(label = "c",
       data  = b1r2,
       color = RUN_COLORS$B1R2,
       y_max =  600, y_step = 100,
       box = c(
         "B1 Run 2 stable phase",
         "(mean s = \u20130.05)",
         paste0(fmt_count(length(b1r2)), " fixed deleterious"),
         sprintf("%.1f%% in ENZ + SBZ", 100 * sum(b1r2 <= SBZ_MAX) / length(b1r2))
       ),
       show_zone_legend = FALSE),

  list(label = "d",
       data  = b1r3,
       color = RUN_COLORS$B1R3,
       y_max =  700, y_step = 100,
       box = c(
         "B1 Run 3 stable phase",
         "(mean s = \u20130.01)",
         paste0(fmt_count(length(b1r3)), " fixed deleterious"),
         sprintf("%.1f%% in ENZ + SBZ", 100 * sum(b1r3 <= SBZ_MAX) / length(b1r3))
       ),
       show_zone_legend = FALSE)
)


# ---------- Panel-building function ------------------------------------------
build_panel <- function(p) {
  # Restrict to displayed range so hist() doesn't reject out-of-bound values
  s_in_range <- p$data[p$data <= X_MAX]
  h <- hist(s_in_range, breaks = BIN_EDGES, plot = FALSE)
  bars <- data.frame(
    xmin = h$breaks[-length(h$breaks)],
    xmax = h$breaks[-1],
    y    = h$counts
  )

  # Zones as a 2-row data frame so legend keys can be generated in panel (a).
  # We use ASCII keys for the fill scale (some R installs choke on unicode in
  # scale keys), and supply pretty labels via the `labels =` argument.
  enz_key <- "ENZ"
  sbz_key <- "SBZ"
  zones <- data.frame(
    xmin = c(0,        ENZ_MAX),
    xmax = c(ENZ_MAX,  SBZ_MAX),
    ymin = c(0,        0),
    ymax = c(Inf,      Inf),
    fill_label = c(enz_key, sbz_key)
  )

  g <- ggplot() +
    # Zone shading
    geom_rect(data = zones,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = fill_label),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(
      name   = NULL,
      values = c(ENZ = COL_ENZ_FILL, SBZ = COL_SBZ_FILL),
      breaks = c("ENZ", "SBZ"),
      labels = c(
        "Effectively neutral zone (|s| < 0.0005)",
        "Selection breakdown zone (0.0005 \u2264 |s| \u2264 0.0025)"
      ),
      guide  = if (isTRUE(p$show_zone_legend))
                 guide_legend(override.aes = list(alpha = 0.5)) else "none"
    ) +

    # Histogram bars (per-run color)
    geom_rect(data = bars,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y),
              fill = p$color, color = "black", linewidth = 0.15,
              alpha = 0.85) +

    # Zone boundary lines
    geom_vline(xintercept = ENZ_MAX, linetype = "dashed",
               color = COL_ENZ_LINE, linewidth = 0.4, alpha = 0.8) +
    geom_vline(xintercept = SBZ_MAX, linetype = "dashed",
               color = COL_SBZ_LINE, linewidth = 0.4, alpha = 0.8) +

    # Info box — top-right, anchored well inside the white space beyond the
    # SBZ shading (|s| > 0.0025) so it doesn't overlap the colored zones.
    # Right edge anchored at 0.998 so it lines up with the zone-legend
    # bounding box below it in panel (a).
    annotate("label",
             x = X_MAX * 0.998, y = p$y_max * 0.97,
             label = paste(p$box, collapse = "\n"),
             hjust = 1, vjust = 1,
             size = 11 / .pt,
             family = JOMB_FONT,
             label.r = unit(0.12, "lines"),
             label.padding = unit(0.38, "lines"),
             fill = "white", color = "black",
             label.size = 0.3,
             lineheight = 1.1) +

    scale_x_continuous(
      name   = expression(bold("|s| (selection coefficient magnitude)")),
      limits = c(0, X_MAX),
      breaks = seq(0, X_MAX, by = 0.001),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
      name   = expression(bold("Fixed deleterious mutations")),
      limits = c(0, p$y_max),
      breaks = seq(0, p$y_max, by = p$y_step),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    theme_jomb()

  # Zone legend only in panel (a). Implemented as two manual annotation
  # blocks (colored swatch + text) placed in the upper-right white space
  # below the info box. This gives us pixel-perfect control over placement,
  # which the auto-positioned ggplot legend can't guarantee when the
  # available white space is narrow.
  if (isTRUE(p$show_zone_legend)) {
    g <- g + theme(legend.position = "none")

    # Legend lives in the rectangle: x in [0.00275, 0.00498], y in [0.46, 0.62]*y_max
    leg_x_left   <- X_MAX * 0.555
    leg_x_swatch <- X_MAX * 0.575
    leg_x_text   <- X_MAX * 0.605
    leg_x_right  <- X_MAX * 0.998

    leg_y_top    <- p$y_max * 0.625
    leg_y_bot    <- p$y_max * 0.420
    leg_y_row1   <- p$y_max * 0.560
    leg_y_row2   <- p$y_max * 0.475
    swatch_h     <- p$y_max * 0.035

    g <- g +
      # Surrounding box
      annotate("rect",
               xmin = leg_x_left, xmax = leg_x_right,
               ymin = leg_y_bot,  ymax = leg_y_top,
               fill = "white", color = "black", linewidth = 0.35) +
      # ENZ swatch + label
      annotate("rect",
               xmin = leg_x_swatch - X_MAX * 0.014,
               xmax = leg_x_swatch + X_MAX * 0.012,
               ymin = leg_y_row1 - swatch_h, ymax = leg_y_row1 + swatch_h,
               fill = COL_ENZ_FILL, color = "black", linewidth = 0.2,
               alpha = 0.7) +
      annotate("text",
               x = leg_x_text, y = leg_y_row1,
               label = "ENZ (|s| < 0.0005)",
               hjust = 0, vjust = 0.5,
               size = (JOMB_BASE_SIZE - 4) / .pt,
               family = JOMB_FONT) +
      # SBZ swatch + label
      annotate("rect",
               xmin = leg_x_swatch - X_MAX * 0.014,
               xmax = leg_x_swatch + X_MAX * 0.012,
               ymin = leg_y_row2 - swatch_h, ymax = leg_y_row2 + swatch_h,
               fill = COL_SBZ_FILL, color = "black", linewidth = 0.2,
               alpha = 0.7) +
      annotate("text",
               x = leg_x_text, y = leg_y_row2,
               label = "SBZ (0.0005 \u2264 |s| \u2264 0.0025)",
               hjust = 0, vjust = 0.5,
               size = (JOMB_BASE_SIZE - 4) / .pt,
               family = JOMB_FONT)
  } else {
    g <- g + theme(legend.position = "none")
  }

  g
}


# ---------- Assemble figure --------------------------------------------------
plots <- lapply(panels, build_panel)

# Use patchwork's auto-tagging ("a", "b", "c", "d"). Each tag is positioned
# in the top margin of its panel, well above the y-axis title (which lives
# in the left axis area). plot.tag.position = c(0, 1) places the tag at the
# very upper-left of the subplot's device area; the generous top margin
# below leaves room so it doesn't collide with the y-axis title.
# ---------- Assemble figure --------------------------------------------------
plots <- lapply(panels, build_panel)

# cowplot::plot_grid places each subplot's label (a/b/c/d) in the white
# space surrounding its upper-left corner — outside the plot area, in the
# margin to the left of the y-axis title and above the top tick. This is
# the placement style used in your fitness-trajectory figure.
fig <- plot_grid(
  plotlist        = plots,
  ncol            = 2,
  labels          = c("a", "b", "c", "d"),
  label_size      = JOMB_BASE_SIZE + 5,
  label_fontface  = "bold",
  label_fontfamily = JOMB_FONT,
  # label_x = 0 anchors the label at the very left edge of each subplot's
  # device area, in the white-space column to the left of the rotated
  # y-axis title. label_y = 1 puts the top of the label at the top edge.
  label_x         = 0,
  label_y         = 1,
  hjust           = -0.2,
  vjust           = 1.4,
  align           = "hv",
  axis            = "tblr"
)


# ---------- Save -------------------------------------------------------------
# 14 x 9 in at 300 dpi, white background — sized for a 2x2 grid where each
# panel has roughly the proportions of the v14 single-panel plots.
ggsave(OUTPUT_PNG, plot = fig,
       width = 14, height = 9, units = "in", dpi = 300, bg = "white")
ggsave(OUTPUT_PDF, plot = fig,
       width = 14, height = 9, units = "in",            bg = "white")

message("Saved: ", OUTPUT_PNG)
message("Saved: ", OUTPUT_PDF)

# =============================================================================
# STYLE NOTES (matches v14)
# =============================================================================
#   Font          : Arial; axis titles bold 15 pt, ticks 14 pt
#   Axes          : black 0.6 pt lines, 4 pt ticks, no gridlines
#   Panel tags    : bold Arial 18 pt, outside upper-left of each panel
#   Per-panel hue : red / blue / green / purple (matches comparison plot)
#   Zone shading  : ENZ pale green #C8E6B8, SBZ pale orange #FFCC99
#   Zone borders  : dashed dark green / dark orange at zone boundaries
#   Output        : 14 x 9 in at 300 dpi, white background
# =============================================================================
