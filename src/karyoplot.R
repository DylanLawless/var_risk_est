# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("karyoploteR",force = TRUE)

# genome wide IEI -----
library(karyoploteR)
library(biomaRt)
library(readr)
library(dplyr)
library(viridisLite)
library(scales)

# Read input data
df <- read_tsv("../output/VarRiskEst_PanelAppRex_ID_398_gene_variants.tsv", show_col_types = FALSE)

# Get unique gene-inheritance pairs
genes_df <- df %>%
  dplyr::select(genename, Inheritance) %>%
  distinct()

# Colour map ----
# colour_map <- c("AD" = "#ff0000", "AR" = "#ff7400", "XL" = "#ffc100")
colour_map <- c("AD" = "#5a00e1", "AR" = "#ff7400", "XL" = "#ffc100")
genes_df$colour <- colour_map[genes_df$Inheritance]

# Get gene coordinates via biomaRt
# mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",  mirror = "useast")
# mart 
# 
# coords <- getBM(
#   attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
#   filters = "hgnc_symbol",
#   values = genes_df$genename,
#   mart = mart
# )
# 
# saveRDS(coords, file="../output/VarRiskEst_PanelAppRex_ID_398_martcoords.Rds")
coords <- readRDS( file="../output/VarRiskEst_PanelAppRex_ID_398_martcoords.Rds")

# Filter valid chromosomes only
coords <- coords %>%
  dplyr::filter(chromosome_name %in% c(1:22, "X", "Y")) %>%
  mutate(chr = paste0("chr", chromosome_name))

# Merge inheritance and colour info
gene_locs <- genes_df %>%
  left_join(coords, by = c("genename" = "hgnc_symbol")) %>%
  dplyr::filter(!is.na(start_position))

# Create modified plot parameters
pp <- getDefaultPlotParams(plot.type = 1)
pp$ideogramheight <- 50   # thicker chromosomes
pp$data1height <- 20      # room for ticks
pp$margin <- 0.01
pp$ideogramlateralmargin <- 0.1
pp$bottommargin <- 50
pp$topmargin <- 50

# Plot karyotype ----

pdf("../images/locusplot_karyoplt.pdf", width = 12, height = 4.2)
kp <- plotKaryotype(genome = "hg38", plot.type = 1, plot.params = pp)

# Plot gene ticks
for (i in seq_len(nrow(gene_locs))) {
  mid <- (gene_locs$start_position[i] + gene_locs$end_position[i]) / 2
  kpSegments(kp,
             chr = gene_locs$chr[i],
             x0 = mid, x1 = mid,
             y0 = -4.5, # end hieght
             y1 = .5, # start height
             col = gene_locs$colour[i],
             lwd = 2)
}

# Add colour legend (overlayed in bottom right)
legend("bottomright",
       legend = c("Autosomal Dominant", "Autosomal Recessive", "X-linked"),
       # col = c("#ff0000", "#ff7400", "#ffc100"),
       col = c("#5a00e1", "#ff7400", "#ffc100"),
       lwd = 6,
       bg = "white",
       box.col = NA,
       cex = 1,
       inset = 0.02)
dev.off()

