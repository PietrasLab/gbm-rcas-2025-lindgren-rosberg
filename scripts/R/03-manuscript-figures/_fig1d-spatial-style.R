# Shared Figure-1D spatial plotting style (from .spatialFeaturePlot in
# scripts/R/00-0-source-functions.R). Source this in any 30-* spatial script that
# uses a gradient colour scale, so all maps share the canonical look:
#   palette : grayish zero -> pale turquoise -> yellow -> orange -> red -> pink -> purple
#   bg      : #F5F3EE (paler than the #E4E3DB tissue/zero), on panel + plot + strip
# Masks/outlines should be drawn ON TOP with a slightly wider linewidth.
suppressPackageStartupMessages(library(ggplot2))
FIG1D_PAL <- c("#E4E3DB", "#BFE8DA", "#FEE08B", "#F7B36B", "#EA8969", "#E07BA3", "#9E3D8D")
FIG1D_BG  <- "#F5F3EE"
MASK_LWD  <- 0.55   # outline linewidth (wider than the old 0.2-0.25)

# thin, short colourbar so continuous legends stay small (independent of legend.key.size,
# which we shrink for DISCRETE keys below)
.fig1d_cbar <- guide_colourbar(barheight = grid::unit(0.3, "lines"), barwidth = grid::unit(5, "lines"),
                               title.position = "top", ticks.colour = NA)
fig1d_color <- function(limits = NULL, name = NULL)
  scale_color_gradientn(colours = FIG1D_PAL, limits = limits, oob = scales::squish, guide = .fig1d_cbar,
                        breaks = if (!is.null(limits)) pretty(limits, 4) else waiver(), name = name)
fig1d_fill <- function(limits = NULL, name = NULL)
  scale_fill_gradientn(colours = FIG1D_PAL, limits = limits, oob = scales::squish, guide = .fig1d_cbar,
                       breaks = if (!is.null(limits)) pretty(limits, 4) else waiver(), name = name)
theme_fig1d <- function(base_size = 10)
  theme_void(base_size = base_size) +
  theme(aspect.ratio = 1, legend.position = "bottom",
        # small legends everywhere (discrete niche/domain keys were dwarfing the maps)
        legend.key.size   = grid::unit(0.32, "lines"),
        legend.text       = element_text(size = 5.5),
        legend.title      = element_text(size = 6.5),
        panel.background  = element_rect(fill = FIG1D_BG, colour = NA),
        plot.background   = element_rect(fill = FIG1D_BG, colour = NA),
        strip.background  = element_rect(fill = FIG1D_BG, colour = NA),
        plot.title        = element_text(face = "bold", hjust = 0.5))
