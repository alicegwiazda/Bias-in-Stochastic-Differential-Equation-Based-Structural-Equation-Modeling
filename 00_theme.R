# =============================================================================
# 00_theme.R
# Visual style for all chapter 8 figures. Source after 00_helpers.R.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})


# -----------------------------------------------------------------------------
# Palette
#   COL_ETA       : baseline / oracle / true latents — near-black, restrained
#   COL_ETAHAT    : two-step EFA scores            — desaturated dark red
#   COL_PREDICTED : analytical prediction curve   — neutral grey, dashed
#   COL_REFERENCE : reference lines (true A, zero) — light grey
# -----------------------------------------------------------------------------

COL_ETA       <- "#1F2A44"   # deep navy-black
COL_ETAHAT    <- "#A23E2A"   # burnt brick red
COL_PREDICTED <- "#6F6F6F"   # mid grey
COL_REFERENCE <- "#B8B8B8"   # light grey

THESIS_FILL <- c("eta" = COL_ETA, "etahat" = COL_ETAHAT)
THESIS_LABS <- c("eta" = "True latents",
                 "etahat" = "EFA factor scores",
                 "predicted" = "Analytical prediction")


# -----------------------------------------------------------------------------
# thesis_theme()
#   ggplot2 theme designed to sit cleanly next to a two-column serif body.
#   Serif font, restrained gridlines, no boxes, subdued caption-grey subtitles.
# -----------------------------------------------------------------------------
thesis_theme <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "serif") +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 1,
                                        margin = margin(b = 4)),
      plot.subtitle      = element_text(colour = "gray35", size = base_size - 1,
                                        margin = margin(b = 8)),
      plot.caption       = element_text(colour = "gray45", size = base_size - 2,
                                        hjust = 0, margin = margin(t = 6)),
      axis.title         = element_text(size = base_size),
      axis.text          = element_text(colour = "gray25"),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "gray92", linewidth = 0.3),
      legend.position    = "top",
      legend.title       = element_blank(),
      legend.text        = element_text(size = base_size - 1),
      legend.key.width   = unit(14, "pt"),
      legend.margin      = margin(b = 2),
      strip.text         = element_text(face = "bold", size = base_size),
      plot.margin        = margin(6, 8, 6, 8)
    )
}


# -----------------------------------------------------------------------------
# save_thesis_plot()
#   Standardised export. Default sizing assumes a two-column ACM-style layout:
#   single-column width ~= 3.4 in, double-column ~= 7.0 in.
# -----------------------------------------------------------------------------
save_thesis_plot <- function(plot, filename,
                             width  = 3.4,
                             height = 2.6,
                             dpi    = 600) {
  ggsave(filename, plot = plot,
         width = width, height = height,
         units = "in", dpi = dpi)
  invisible(filename)
}