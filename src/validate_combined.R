library(patchwork)

# get validation result plots
# Be careful here as overlapping environment variables could get confused 
# Likewise, note that variables would be resued so see individual scripts for datasets

# dominant example
plot_nfkb1 <- local({
  source("validate_nfkb1.R")
  p_nfkb1_bayes
})

# rm(list = setdiff(ls(), "plot_nfkb1"))

# recessive example
plot_cftr <- local({
  source("validate_cftr.R")
  p_cftr_bayes
})

# Now, if you want to clear other variables:
# rm(list = setdiff(ls(), c("plot_nfkb1", "plot_cftr")))

# plot_nfkb1
# plot_cftr

patch1 <- patch1  + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
patch1

ggsave("../images/validation_studies_bayesian_adjusted_estimates.png", plot = patch1, width = 9, height = 6)
