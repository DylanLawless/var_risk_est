# HPC large import ----
library(readr)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2); theme_set(theme_bw())

# save data ----
df <- readRDS(file = "../data/panel_all_genes_df.Rds")
df_tally <- readRDS(file = "../data/panel_all_genes_df_tally.Rds")

var_count <- df |> nrow()
gene_count <- df$genename |> unique() |> length()
panel_name <- df$name |> unique()
panel_ID <- df$panel_id |> unique()

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


# table ----
library(dplyr)
library(knitr)
library(kableExtra)

df_subset <- df %>%
  select(genename, panel_id, clinvar_clnsig, `pos(1-based)`, HGVSc_VEP, HGVSp_VEP, Inheritance, occurrence_prob) %>% arrange(genename, `pos(1-based)`)

df_subset$clinvar_clnsig <- gsub("_", " ", df_subset$clinvar_clnsig)

df_subset_print <- df_subset %>%
  select(genename, panel_id, clinvar_clnsig, `pos(1-based)`, HGVSc_VEP, HGVSp_VEP, Inheritance, occurrence_prob) %>%
  # Rename columns for publication-quality output
  rename(
    "Gene" = genename,
    "Panel ID" = panel_id,
    "ClinVar Clinical Significance" = clinvar_clnsig,
    "GRCh38 Pos" = `pos(1-based)`,
    "HGVSc (VEP)" = HGVSc_VEP,
    "HGVSp (VEP)" = HGVSp_VEP,
    "Inheritance" = Inheritance,
    "Occurrence Probability" = occurrence_prob
  ) %>%
  # Format numeric column to retain precision (9 decimal places)
  mutate(`Occurrence Probability` = sprintf("%.9f", `Occurrence Probability`))

library(knitr)
library(kableExtra)

kable(head(df_subset_print), 
      format = "latex", 
      booktabs = TRUE,
      escape = FALSE,
      caption = paste0("Example subset of the main results for ", gene_count, " genes of PanelAppRex's panel: (ID ", panel_ID, ") ", panel_name, ". ", "``ClinVar Significance'' indicates the pathogenicity classification assigned by ClinVar, while ``Occurrence Prob'' represents our calculated probability of observing the corresponding variant class for a given phenotype.")) %>%
  kable_styling(latex_options = c("scale_down"))

# Summary of PAR -----
source("panelapprex_import.R")
df_par <- df_par |> filter(panel_id == panel_id ) |> slice(1)
df_str <- paste(capture.output(print(df_par)), collapse = "\n")
cat(df_str)


# Result
paste0("We report the probability of disease observation for a total of ",
       var_count, " ClinVar variant classifications in ", gene_count, " genes.")

# Discussion
paste0("In this study we focused on reporting the probability of disease observation for the  disease gene panel (ID ", panel_ID, ") ", panel_name, ". We report a total of ",
var_count, " ClinVar variant classifications in ", gene_count, " genes. This demonstrates in a thoruough complete example how we can additionally run through all gene-disease combinations genome wide.")





