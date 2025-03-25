
# Import PanelAppRex Rds format ---
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_par <- readRDS(file= path_PanelAppData_genes_combined_Rds)

# Set Inheritance, if unknown default to AR for safety ----
df_par <- df_par %>%
  mutate(Inheritance = case_when(
    grepl("X[- ]?LINKED", mode_of_inheritance, ignore.case = TRUE) ~ "XL",
    grepl("both", mode_of_inheritance, ignore.case = TRUE) ~ "AR",
    grepl("biallelic", mode_of_inheritance, ignore.case = TRUE) ~ "AR",
    grepl("monoallelic", mode_of_inheritance, ignore.case = TRUE) ~ "AD",
    TRUE ~ "AR"
  ))

