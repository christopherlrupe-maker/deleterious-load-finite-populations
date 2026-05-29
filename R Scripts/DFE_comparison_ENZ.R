# =====================================================================
# DFE Comparison Plot (ENZ Shaded) — v14 JoMB Styling
# =====================================================================
# Shows the gamma distributions of fitness effects for H&C's baseline
# and three B1 runs isolating mean fitness effects. The effectively
# neutral zone (ENZ, |s| < 0.0005) is shaded in pale green.
#
# Style matches the v14 JoMB plot generator family:
#   - Arial font, bold axis titles at 15 pt
#   - No gridlines, axis lines on left and bottom only
#   - Top-left legend, no border
#   - White background, 300 dpi output
#   - 4-color palette preserved for the 4 run curves
#   - ENZ shading and dashed boundary line preserved
#
# Title and subtitle are removed (descriptive text belongs in the
# manuscript caption, not on the figure itself).
#
# Requires: ggplot2
#   install.packages("ggplot2")
# =====================================================================

library(ggplot2)

# ---------- JoMB style constants (matching v14) ----------------------
JOMB_FONT      <- "Arial"
JOMB_BASE_SIZE <- 15

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

      # Legend: top-left inside panel, no border
      legend.position      = c(0.03, 0.97),
      legend.justification = c("left", "top"),
      legend.background    = element_blank(),
      legend.key           = element_rect(fill = "white", color = NA),
      legend.title         = element_blank(),
      legend.text          = element_text(size = base_size - 1),
      legend.spacing.y     = unit(2, "pt"),
      legend.margin        = margin(2, 4, 2, 4),

      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin      = margin(10, 14, 6, 10)
    )
}

# ---------- Parameters ----------------------------------------------
# Each run uses shape = 0.5 (H&C's value). Gamma DFE is parameterized
# by shape (a) and scale (theta); mean = shape * scale, so for a given
# absolute mean fitness effect, scale = |mean_s| / shape.
shape <- 0.5

runs <- data.frame(
  label  = c("H&C baseline (mean s = -0.5)",
             "B1 Run 1 (mean s = -0.1)",
             "B1 Run 2 (mean s = -0.05)",
             "B1 Run 3 (mean s = -0.01)"),
  mean_s = c(0.5, 0.1, 0.05, 0.01),
  color  = c("#C9302C",   # red
             "#1F4E79",   # blue
             "#1F6F2F",   # green
             "#5E2D8A"),  # purple
  stringsAsFactors = FALSE
)
runs$scale <- runs$mean_s / shape

# ---------- Curve data --------------------------------------------------
# x range covers the deleterious side of the DFE down to s = -0.010.
x_seq <- seq(-0.010, 0, length.out = 2000)

curve_data <- do.call(rbind, lapply(seq_len(nrow(runs)), function(i) {
  density <- dgamma(abs(x_seq), shape = shape, scale = runs$scale[i])
  data.frame(
    s       = x_seq,
    density = density,
    run     = factor(runs$label[i], levels = runs$label),
    stringsAsFactors = FALSE
  )
}))

# ---------- Zone boundary --------------------------------------------
ENZ_THRESHOLD <- 0.0005   # |s| < 1/(2N) for N=1000

run_colors <- setNames(runs$color, runs$label)

# ---------- Build the plot --------------------------------------------
p1 <- ggplot(curve_data, aes(x = s, y = density, color = run)) +

  # Shade ENZ in pale green and add it to the legend via aes(fill = ...).
  # We build a one-row data frame with the zone bounds and map the fill to a
  # label string; scale_fill_manual below assigns the actual pale green color.
  geom_rect(data = data.frame(
              xmin = -ENZ_THRESHOLD, xmax = 0,
              ymin = 0, ymax = Inf,
              fill_label = "Effectively neutral zone (|s| < 0.0005)"),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = fill_label),
            inherit.aes = FALSE,
            alpha = 0.5) +
  scale_fill_manual(
    name = NULL,
    values = setNames("#C8E6B8",
                      "Effectively neutral zone (|s| < 0.0005)"),
    guide = guide_legend(order = 2,
                         override.aes = list(alpha = 0.5))
  ) +

  # Dashed boundary line at the LEFT edge of the ENZ (at s = -0.0005)
  geom_vline(xintercept = -ENZ_THRESHOLD,
             linetype = "dashed", color = "#2E7D32",
             linewidth = 0.4, alpha = 0.8) +
  # Solid black line at the RIGHT edge of the ENZ (s = 0)
  geom_vline(xintercept = 0,
             linetype = "solid", color = "black",
             linewidth = 0.6) +

  # Plot the gamma curves — slightly thicker than trajectory lines because
  # these are smooth analytical curves, not noisy simulation traces.
  geom_line(linewidth = 0.7) +

  # Custom color palette for the 4 runs
  scale_color_manual(values = run_colors, name = NULL,
                     guide = guide_legend(order = 1)) +

  # Axes — bold, sans-serif via theme_jomb()
  scale_x_continuous(
    name   = expression(bold("Selection coefficient (s) of deleterious mutations")),
    limits = c(-0.010, 0),
    breaks = seq(-0.010, 0, by = 0.002),
    expand = expansion(mult = c(0, 0.005))   # zero left padding so curves touch y-axis
  ) +
  scale_y_continuous(
    name   = expression(bold("Probability density")),
    limits = c(0, 250),
    expand = expansion(mult = c(0, 0))
  ) +

  # Apply the v14-matching theme (after all scales/aesthetics)
  theme_jomb() +
  # Keep the legend INSIDE the panel, flush against the y-axis (matches the
  # v14 style of other plots). Earlier worry about overlap with the axis was
  # misplaced: at x = 0, the legend sits to the right of where the y-axis
  # line is drawn, not on top of it. legend.position = c(0, 1) with
  # justification c("left", "top") anchors it tight against the upper-left
  # corner of the plotting area.
  theme(legend.position = c(0.03, 1),
        legend.justification = c("left", "top")) +

  # Annotation for ENZ — bold green text at the top of the shaded zone
  annotate("text",
           x = -ENZ_THRESHOLD / 2,
           y = 8,                   # bottom of the green-shaded zone
           label = "ENZ",
           color = "#2E7D32",
           size  = 5,
           fontface = "bold",
           family = JOMB_FONT)

# ---------- Save ---------------------------------------------------------
# 10 x 6 inches at 300 dpi, white background — slightly wider than the
# v14 trajectory output (8 x 6) because this plot has a long y-axis label
# legend that benefits from extra horizontal room.
ggsave("DFE_comparison_ENZ_only.png",
       plot   = p1,
       width  = 10,
       height = 6,
       dpi    = 300,
       bg     = "white")

print(p1)
cat("Saved: DFE_comparison_ENZ_only.png\n")

# =====================================================================
# STYLE NOTES (matches v14)
# =====================================================================
#   Font         : Arial; axis titles bold 15 pt, ticks/legend 14 pt
#   Axes         : black 0.6 pt lines, 4 pt ticks, no gridlines
#   Legend       : top-left (c(0.03, 0.97)), no border
#   Output       : 10 x 6 in at 300 dpi, white background
#   ENZ shading  : pale green #C8E6B8, dashed boundary at -0.0005
#   Curve colors : red / blue / green / purple — 4 distinct runs
#   Line width   : 0.7 pt (smooth analytical curves)
# =====================================================================
