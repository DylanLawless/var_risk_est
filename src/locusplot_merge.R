source("locusplot_nfkb1.R")
source("locusplot_cftr.R")
source("karyoplot.R")

library(patchwork)
library(magick)
library(grid)
library(ggplot2)
library(pdftools)

# Read in PDFs as high-res images
img0 <- image_read_pdf("../images/locusplot_karyoplt.pdf", density = 300)
img1 <- image_read_pdf("../images/locusplot_nfkb1.pdf", density = 300)
img2 <- image_read_pdf("../images/locusplot_cftr.pdf", density = 300)

# Convert to grobs
g0 <- rasterGrob(img0, width = unit(1, "npc"), height = unit(1, "npc"), just = "centre")
g1 <- rasterGrob(img1, width = unit(1, "npc"), height = unit(1, "npc"), just = "centre")
g2 <- rasterGrob(img2, width = unit(1, "npc"), height = unit(1, "npc"), just = "centre")

design <- "
00
12"

# Wrap and layout
p_combined <- wrap_elements(g0) / (wrap_elements(g1) + wrap_elements(g2)) +
  # plot_layout(design = design, guides = "collect") +
  plot_layout(nrow = 2, heights = c(2, 2.5), guides = "collect") +
  plot_annotation(tag_levels = "A")

p_combined

# Save output (PNG or PDF)
ggsave("../images/karyo_locusplot_merged.pdf", p_combined, width = 8, height = 6, device = cairo_pdf)

 # # locus only -----
# design <- "
# 12"
# 
# # Wrap and layout
# p_combined <- wrap_elements(g1) + wrap_elements(g2) +
#   plot_layout(design = design, guides = "collect") +
#   plot_annotation(tag_levels = "A")
# 
# p_combined
# 
# # Save output (PNG or PDF)
# ggsave("../images/locusplot_merged.pdf", p_combined, width = 6, height = 12, device = cairo_pdf)
