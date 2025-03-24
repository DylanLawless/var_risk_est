# HPC large import ----
library(readr)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2); theme_set(theme_bw())

# save data ----
df <- readRDS(file = "../data/panel_all_genes_df.Rds")
df_tally <- readRDS(file = "../data/panel_all_genes_df_tally.Rds")

gene_count <- df$genename |> unique() |> length()

# Bar plot: Count of ClinVar Clinical Significance (across all genes)
p_count <- df |>
  filter(synth_flag == "FALSE") |>
 ggplot(aes(
  x = stringr::str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  fill = clinvar_clnsig)) +
  geom_bar(color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = paste("Count of ClinVar Clinical Significance (All Genes = ", gene_count, ")")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave("../images/genome_all_genes_clinvar_count.png", plot = p_count, width = 10, height = 4)
