library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Construct a toy dataset
toy_variants <- tibble::tribble(
  ~chr,   ~pos,    ~ref, ~alt, ~genotype, ~tested, ~in_ref,
  "chr1", 101,     "A",  "G",  0,         TRUE,    TRUE,   # known, tested, hom-ref
  "chr1", 102,     "C",  "T",  1,         TRUE,    TRUE,   # known, tested, het
  "chr1", 103,     "G",  "A",  2,         TRUE,    TRUE,   # known, tested, hom-alt
  "chr1", 104,     "T",  "C",  NA,        FALSE,   TRUE,   # known, untested
  "chr2", 201,     "A",  "T",  1,         TRUE,    FALSE,  # novel, observed het
  "chr2", 202,     "G",  "C",  2,         TRUE,    FALSE,  # novel, observed hom-alt
  "chr2", 204,     "T",  "G",  0,         TRUE,    FALSE   # novel position, observed as reference
)


# Classification column
toy_variants <- toy_variants %>%
  mutate(
    classification = case_when(
      in_ref & tested & !is.na(genotype) ~ "observed_known",
      in_ref & !tested                   ~ "untested_known",
      !in_ref & tested & !is.na(genotype) ~ "novel_observed",
      TRUE                               ~ "other"
    ),
    genotype_label = factor(genotype, levels = c(0, 1, 2), labels = c("hom_ref", "het", "hom_alt"))
  )

# View dataset
print(toy_variants)

# Add 'missing' as a level before using replace_na
summary_counts <- toy_variants %>%
  count(classification, genotype_label, name = "n") %>%
  mutate(genotype_label = fct_expand(genotype_label, "missing")) %>%
  replace_na(list(genotype_label = "missing"))

print(summary_counts)

# Define colour scale for genotypes including 'missing'
geno_colours <- c(
  "hom_ref"  = "orange",  # light grey for 0 dosage
  "het"      = "red1",  # blue for 1 dosage
  "hom_alt"  = "red4",  # darker blue for 2 dosage
  "missing"  = "navy"   # very light grey for unsequenced/missing
)

# Plot
ggplot(summary_counts, aes(x = classification, y = n, fill = genotype_label)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = geno_colours) +
  labs(title = "Variant classification observation and genotype dosage",
       x = "Classification",
       y = "Count",
       fill = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# make VCF ----

# Function to convert numeric genotype to VCF GT string
geno_to_gt <- function(gt) {
  if (is.na(gt)) return("./.")  # uncalled
  if (gt == 0) return("0/0")
  if (gt == 1) return("0/1")
  if (gt == 2) return("1/1")
  return("./.")
}

# Construct VCF body
vcf_body <- toy_variants %>%
  filter(tested == TRUE) %>%  # exclude untested lines
  mutate(
    ID = ".",
    QUAL = ".",
    FILTER = "PASS",
    INFO = ".",
    FORMAT = "GT",
    sample = sapply(genotype, geno_to_gt)
  ) %>%
  select(`#CHROM` = chr, POS = pos, ID, REF = ref, ALT = alt,
         QUAL, FILTER, INFO, FORMAT, `patient1` = sample)

# VCF header lines
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##source=toyVCFGenerator",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpatient1"
)

# Write VCF to file
vcf_file <- "../data/toy_patient.vcf"
readr::write_lines(vcf_header, vcf_file)
readr::write_tsv(vcf_body, vcf_file, append = TRUE, col_names = FALSE)

# make checklist ----
checklist <- tibble::tribble(
  ~chr,   ~pos, ~ref, ~alt,
  "chr1", 101,  "A",  "G",
  "chr1", 102,  "C",  "T",
  "chr1", 103,  "G",  "A",
  "chr1", 104,  "T",  "C"  # this one is missing in the VCF (simulates untested)
)

readr::write_tsv(checklist, "../data/checklist.tsv", col_names = FALSE)
