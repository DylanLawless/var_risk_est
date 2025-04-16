

# our data ----

library(ggplot2);theme_set(theme_bw())
library(dplyr)

PANEL <- 398

# PanelAppRex ----
# Rds format
# path_data <- "~/web/PanelAppRex/data"
# path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
# df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
# colnames(df_core)[colnames(df_core) == 'entity_name'] <- 'Genetic defect'
# df_core <- df_core |> filter(panel_id == PANEL)

# VarRiskEst ----
# source("panelapprex_import.R")

# Load gene tally and variant data for the specified panel using the PANEL variable
varRisEst_gene <- readRDS(file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_tally.Rds"))
varRisEst_var  <- readRDS(file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_variants.Rds"))

df <- varRisEst_var |> dplyr::filter(genename == "RAG1") |> dplyr::filter(clinvar_clnsig == "Pathogenic")
df$CHR <- 11
df$BP <- df$`pos(1-based)`

bp_min <- df$BP |> min()
bp_max <- df$BP |> max()

# variants <- df %>%
#   mutate(variant = paste0("chr", CHROM, ":", POS)) %>%
#   pull(variant)


# https://pan.ukbb.broadinstitute.org/blog/2020/10/29/ld-release
# 
# https://github.com/bulik/ldsc

library(ggplot2)


# Instead of using a system pipe with head (which may fail on some systems),
# we can read directly from the gzipped file using gzfile() and the nrows parameter.
# ldscore <- read.table(gzfile("~/Desktop/UKBB.ALL.ldscore/UKBB.EUR.l2.ldscore.gz"),
                      # header = TRUE, stringsAsFactors = FALSE, nrows = 10000)

ldscore <- read.table(gzfile("~/Desktop/UKBB.ALL.ldscore/UKBB.EUR.l2.ldscore.gz"),
                      header = TRUE, stringsAsFactors = FALSE)
ldscore$Index <- 1:nrow(ldscore)

ldscore <- ldscore |> dplyr::filter(CHR == 11)
ldscore$Index <- 1:nrow(ldscore)

ldscore <- ldscore |> dplyr::filter(BP > 36074287)
ldscore <- ldscore |> dplyr::filter(BP < 37074287)

# ldscore <- ldscore |> dplyr::filter(BP >= bp_min/2)
# ldscore <- ldscore |> dplyr::filter(BP < bp_max*2)

ggplot(ldscore, aes(x = BP, y = L2)) +
  geom_line() +
  labs(x = "Variant Index", y = "L2 Score", title = "UKBB LD Scores (First 1000 Variants)")



# LDlinkR ----

# install.packages("LDlinkR")
library(LDlinkR)
# browseVignettes("LDlinkR")


# Set your LDlink API token here
api_token <- "***"

# Example: Get pairwise r2 between two variants
# snp1 <- "chr1:1692321:T:C"
# snp2 <- "1:1695574:G:T"
# snp1 <- "chr1:1692321"
# snp2 <- "chr1:1695574"
# 
# pair_ld <- LDpair(var1 = snp1, var2 = snp2, pop = "EUR", token = api_token)
# print(pair_ld)


library(LDlinkR)

# Your list of variants (adjust alleles if necessary)
# variants <- c("chr1:752721", "chr1:754182", "chr1:768448", "chr1:838555")

# Query pairwise LD (R2) for these variants in European (CEU) population, GRCh37 build
ld_matrix <- LDmatrix(
  snps = variants,
  pop = "CEU",
  r2d = "r2",
  token = api_token,  # or specify your token as a string
  genome_build = "grch38"
)

print(ld_matrix)
