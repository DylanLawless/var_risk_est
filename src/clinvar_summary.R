# see the iei genetics code

# clinvar summary ----

# library(tidyverse)
# library(readr)
# library(ggplot2)
library(stringr)

# Assume 'df' is already loaded. Rename clinvar_clnsig to clinvar_clnsig.
# df <- df %>% rename(clinvar_clnsig = clinvar_clnsig)

# Use the same score_map as before.
score_map <- c(
			     "Pathogenic" = 5,
				   "Pathogenic/Likely_pathogenic" = 4,
				   "not_provided" = 0,
				     "Uncertain_significance" = 0,
				     "Conflicting_interpretations_of_pathogenicity" = 2,
					   "Likely_pathogenic" = 4,
					   "Likely_benign" = -4,
					     "Pathogenic|drug_response" = 5,
					     "drug_response" = 0,
						   "Benign/Likely_benign" = -5
						 )

# Function to compute rank score for a clinvar_clnsig string.
rank_clin_sig <- function(clin_str) {
	  unified <- str_replace_all(clin_str, ";", "/")
  terms <- str_split(unified, "/")[[1]] %>% str_trim()
    scores <- sapply(terms, function(term) {
						     key <- names(score_map)[str_to_lower(names(score_map)) == str_to_lower(term)]
							     if (length(key) > 0) {
									       score_map[[key]]
							     } else {
									       0
								     }
							   })
    mean(scores)
}

# Compute counts per clinvar_clnsig.
final_results <- df %>%
	  group_by(clinvar_clnsig) %>%
	    summarise(total = n(), .groups = "drop")

	# Compute rank for each unique clinvar_clnsig.
	ranked_terms <- final_results %>%
		  mutate(Rank = sapply(clinvar_clnsig, rank_clin_sig))

	  # For demonstration, use the clinical classification as the grouping variable.
	  # p1: plot average rank vs. count of variants per clinvar_clnsig.
	  p1 <- ggplot(ranked_terms, aes(x = Rank, y = total, fill = Rank)) +
		    geom_point(shape = 21, size = 3, color = "black") +
			  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
			    labs(
					     x = "Average Rank Score",
						     y = "Count of Variants",
						     fill = "Avg Rank"
							   ) +
  theme_bw()

  print(p1)


