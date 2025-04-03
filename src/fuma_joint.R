library(tidyverse)
library(patchwork)
library(stringr)

# 1) Import both cluster_2 and cluster_4 data, add a column to label each
data2 <- read_tsv("../data/fuma/FUMA_gene2func604419_var_risk_est_cluster_2/GS.txt") %>%
  mutate(study = "cluster_2")
data4 <- read_tsv("../data/fuma/FUMA_gene2func604403_var_risk_est_cluster_4/GS.txt") %>%
  mutate(study = "cluster_4")

# 2) Combine both into a single data frame
all_data <- bind_rows(data2, data4)

common_theme <- theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# 3) Function to produce a composite plot and save full GeneSet names to file
plotStudy <- function(data, study_label, category) {
  # Filter and transform data
  data <- data %>%
    filter(Category == category) %>%
    mutate(
      proportion = N_overlap / N_genes,
      enrichment = -log10(adjP),
      GeneSet = str_replace_all(GeneSet, "_", " ")
    ) %>%
    filter(adjP < 0.05) %>%
    arrange(adjP) %>%
    slice_head(n = 25) %>%
    # Convert GeneSet to sentence case (full names) before truncation
    mutate(GeneSet = str_to_sentence(GeneSet))
  
  # Save the full GeneSet names (unique) to a file, including study and category in filename
  full_gene_sets <- data %>% pull(GeneSet) %>% unique() %>% paste(collapse = "\n")
  # Replace spaces in study_label with underscore
  study_label_file <- str_replace_all(study_label, " ", "_")
  filename <- paste0("../data/fuma/fum_", study_label_file, "_", category, "_geneset_summary.txt")
  write_lines(full_gene_sets, filename)
  
  # Now truncate GeneSet labels to width 30 characters
  data <- data %>%
    mutate(GeneSet = str_trunc(GeneSet, width = 30, side = "right", ellipsis = "...")) %>%
    # Retain only the row with the largest N_genes for each GeneSet if duplicates exist
    group_by(GeneSet) %>%
    filter(N_genes == max(N_genes)) %>%
    ungroup() %>%
    # Ensure y axis has unique labels in the desired order
    mutate(GeneSet = factor(GeneSet, levels = rev(unique(GeneSet))))
  
  max_prop <- max(data$proportion, na.rm = TRUE)
  max_enrich <- max(data$enrichment, na.rm = TRUE)
  
  # Set colours based on category
  if (category == "GWAScatalog") {
    fill_left <- "#64abaa"
    fill_right <- "#add6d5"
  } else if (category == "Immunologic_signatures") {
    fill_left <- "#ebb382"
    fill_right <- "#f7ca89"
  } else {
    fill_left <- "#107dab"
    fill_right <- "#ff0000"
  }
  
  p_left <- ggplot(data, aes(x = proportion, y = GeneSet)) +
    geom_col(fill = fill_left, color = "black") +
    scale_x_reverse(limits = c(max_prop, 0), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "proportion") +
    common_theme +
    ggtitle(study_label) +
    xlim(c(1, 0))
  
  p_right <- ggplot(data, aes(x = enrichment, y = GeneSet)) +
    geom_col(fill = fill_right, color = "black") +
    scale_x_continuous(limits = c(0, max_enrich), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "-log10\nadjusted p-value") +
    common_theme +
    xlim(c(0, 13))
  
  p_left + p_right +
    plot_layout(ncol = 2, widths = c(1, 1), guides = "collect", axes = "collect")
}

# 4) Create composite plots for each study and category using the combined data
fuma_gwascat2 <- plotStudy(all_data %>% filter(study == "cluster_2"), "Cluster 2 - GWAScatalog", "GWAScatalog")
fuma_gwascat4 <- plotStudy(all_data %>% filter(study == "cluster_4"), "Cluster 4 - GWAScatalog", "GWAScatalog")
fuma_imm2     <- plotStudy(all_data %>% filter(study == "cluster_2"), "Cluster 2 - Immune signature", "Immunologic_signatures")
fuma_imm4     <- plotStudy(all_data %>% filter(study == "cluster_4"), "Cluster 4 - Immune signature", "Immunologic_signatures")

# 5) Arrange the GWAScatalog plots vertically
fuma_gwascat <- fuma_gwascat2 / fuma_gwascat4 + 
  plot_layout(ncol = 1, guides = "collect", axes = "collect")

# 6) Arrange the Immunologic_signatures plots vertically
fuma_imm <- fuma_imm2 / fuma_imm4 + 
  plot_layout( guides = "collect", axes = "collect")

fuma_gwascat
fuma_imm

ggsave("../images/fuma_gwascat.pdf", fuma_gwascat, width = 6, height = 7, device = cairo_pdf)
ggsave("../images/fuma_imm.pdf", fuma_imm, width = 6, height = 7, device = cairo_pdf)

fuma_merge <- fuma_gwascat | fuma_imm  

ggsave("../images/fuma_merge.pdf", fuma_merge, width = 12, height = 6.5, device = cairo_pdf)
