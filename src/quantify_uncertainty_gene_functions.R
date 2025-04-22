
# simulate cohort AN --------------------------------------------------------
cohort_size <- 200
patient_data <- patient_data %>%
  mutate(
    cohort_missing = map_int(known_var_observed, ~ if(.x==-9L)1L else 0L) +
      rbinom(n(), cohort_size-1L, 0.10),
    cohort_AN      = 2L * (cohort_size - cohort_missing)
  ) %>%
  select(-cohort_missing)


# select and score risk_variants -------------------------------------------
risk_variants <- patient_data %>%
  filter((known_var_observed==1 & score>=0) | known_var_observed==-9) %>%
  mutate(
    flag       = if_else(known_var_observed==-9, "missing", "present"),
    class      = if_else(evidence_score>3, "causal", "other"),
    variant_lab= if_else(HGVSp_VEP==".", HGVSc_VEP, HGVSp_VEP),
    alpha      = round(adj_occurrence_prob * cohort_AN) + prior_w,
    beta       = cohort_AN - round(adj_occurrence_prob * cohort_AN) + 1,
    group      = factor(paste(flag, class),
                        levels=c("present causal","present other",
                                 "missing causal","missing other"))
  )

# posterior sampling -------------------------------------------------------
n_sim   <- 10000
# n_sim   <- 5000
# sim_mat <- sapply(seq_len(nrow(risk_variants)),
#                   function(i) rbeta(n_sim,
#                                     risk_variants$alpha[i],
#                                     risk_variants$beta[i]))
# 
# head(sim_mat)

# faster and better memory
sim_mat <- matrix(
  rbeta(n_sim * nrow(risk_variants), 
        rep(risk_variants$alpha, each = n_sim), 
        rep(risk_variants$beta,  each = n_sim)),
  nrow = n_sim, byrow = FALSE
)
# head(sim_mat)

colnames(sim_mat) <- risk_variants$variant_lab

# normalisations ------------------------------------------------------------
# raw for A & B
share_all  <- sim_mat / rowSums(sim_mat)
# causal‐only for C
causal_ix  <- which(risk_variants$class=="causal")
sim_causal <- share_all[, causal_ix, drop=FALSE]
post_share <- colMeans(sim_causal)

risk_variants <- risk_variants %>%
  mutate(
    posterior_share      = if_else(class=="causal",
                                   post_share[variant_lab], 0),
    prob_causal_damaging = posterior_share
  )

# credible intervals for A -------------------------------------------------
ci_present <- risk_variants %>%
  filter(flag=="present") %>%
  transmute(
    variant_lab,
    lower  = quantile(share_all[,variant_lab],0.025),
    median = quantile(share_all[,variant_lab],0.500),
    upper  = quantile(share_all[,variant_lab],0.975)
  )

plot_df <- risk_variants %>%
  left_join(ci_present, by="variant_lab") %>%
  mutate(
    lower    = if_else(flag=="present", lower, posterior_share),
    median   = if_else(flag=="present", median, posterior_share),
    upper    = if_else(flag=="present", upper, posterior_share),
    var_plot = reorder_within(variant_lab, median, group)
  )

# — plot A: raw priors (all variants) —————————————————————————
# prior_df <- risk_variants %>%
  # mutate(var_plot = reorder_within(variant_lab, occurrence_prob, group))
prior_df <- risk_variants %>%
  mutate(var_plot = fct_reorder(variant_lab, score, .desc = FALSE))

p_prior <- ggplot(prior_df,
                  aes(x = occurrence_prob, y = var_plot, fill = score)) +
  geom_point(shape = 21, colour = "black", size = 3, stroke = 0.3) +
  geom_text(aes(label = score),
            hjust = 0, nudge_x = 0.0005, size = 3) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks = c(-5, 0, 5),
                       labels = c("benign", "unknown", "pathogenic")) +
  labs(x = "prior occurrence probability", y = NULL) +
  scale_x_continuous(limits = c(0, 0.0035),
                     breaks = pretty_breaks(n = 2))

# p_prior

# — plot B: full posterior distributions ——————————————————————————
share_df <- as_tibble(share_all) %>%
  mutate(draw=row_number()) %>%
  pivot_longer(-draw, names_to="variant_lab", values_to="share") %>%
  left_join(risk_variants %>% select(variant_lab, group, score),
            by="variant_lab") %>%
  # mutate(var_plot=reorder_within(variant_lab, share, group))
  mutate(var_plot = fct_reorder(variant_lab, score, .desc = FALSE))

p_dist <- ggplot(share_df,
                 aes(x=share, y=var_plot, fill=score)) +
  geom_density_ridges(scale=1, rel_min_height=0.01,
                      colour="grey30", alpha=0.7) +
  facet_grid(group~., scales="free_y", space="free_y",
             labeller=labeller(group=label_wrap_gen(9))) +
  theme(strip.text.y=element_text(size=7)) +
  scale_fill_gradientn(colours=c("navy","lightblue","red"),
                       breaks=c(-5,0,5),
                       labels=c("benign","unknown","pathogenic")) +
  scale_y_reordered() +
  labs(x="posterior share distribution\n(causal & damaging)", y=NULL) +
  scale_x_continuous(limits   = c(0, 1), n.breaks = 3)

# — plot C: p(causal & damaging) ——————————————————————————————

# risk_variants %>% mutate(var_plot=reorder_within(
# variant_lab, prob_causal_damaging, group)
# ),

p_causal <- ggplot(prior_df,
  aes(x=prob_causal_damaging, y=var_plot, fill=score)
) +
  geom_point(shape=21, colour="black", size=3) +
  facet_grid(group~., scales="free_y", space="free_y",
             labeller=labeller(group=label_wrap_gen(9))) +
  theme(strip.text.y=element_text(size=7)) +
  scale_fill_gradientn(colours=c("navy","lightblue","red"),
                       breaks=c(-5,0,5),
                       labels=c("benign","unknown","pathogenic")) +
  scale_y_reordered() +
  labs(x="p(causal & damaging)", y=NULL) +
  scale_x_continuous(limits   = c(0, 1), n.breaks = 3)
  # scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))

# plot D: causal-only density ridges --------------------------------------
share_df_causal <- share_df %>% filter(variant_lab %in% names(post_share))

p_causal_dist <- ggplot(share_df_causal, aes(x=share,y=var_plot,fill=group)) +
  geom_density_ridges(scale=1, rel_min_height=0.01, colour="grey30", alpha=0.7) +
  facet_grid(group~., scales="free_y", space="free_y",
                      labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_manual(values=c("present causal"="forestgreen","missing causal"="orange")) +
  scale_y_reordered() + labs(x="damaging-only \nposterior p(causal & damaging)",y=NULL) +
  theme(legend.position="none") +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))

# plot E: stacked contribution --------------------------------------------
variant_source_df <- risk_variants %>% filter(class=="causal") %>%
  select(flag,prob_causal_damaging) %>% mutate(prob_plot=prob_causal_damaging+1e-2)
source_totals <- variant_source_df %>% group_by(flag) %>% summarise(p_damaging=sum(prob_causal_damaging)) %>% mutate(prop=p_damaging/sum(p_damaging))

# p_source_stack <- ggplot(variant_source_df, aes(x=flag,y=prob_plot,fill=flag)) +
#   geom_col(width=0.6,colour="black",position="stack") +
#   geom_text(data=source_totals, aes(x=flag,y=p_damaging,label=scales::percent(prop,accuracy=0.1)), vjust=-0.3,size=3.5) +
#   scale_fill_manual(values=c("present"="forestgreen","missing"="orange")) +
#   labs(x=NULL,y="total p(causal & damaging)") + theme_bw()

# stacked contribution including zero‑slices for non‑causals
eps <- 1e-2
variant_source_df <- risk_variants %>%
  # include all variants, even non‑causal
  select(flag, score, prob_causal_damaging, variant_lab) %>%
  mutate(
    # add epsilon so zero‑prob variants still get a sliver
    prob_plot = prob_causal_damaging + eps
  )

source_totals <- variant_source_df %>%
  group_by(flag) %>%
  summarise(p_damaging = sum(prob_causal_damaging), .groups = "drop") %>%
  mutate(prop = p_damaging / sum(p_damaging))


variant_source_df <- variant_source_df %>%
  arrange(flag, score) %>%
  group_by(flag) %>%
  mutate(
    ymax = cumsum(prob_plot),
    ymin = ymax - prob_plot,
    ymid = (ymin + ymax) / 2
  ) %>%
  ungroup()

p_source_stack <- ggplot(variant_source_df,
                         aes(x = flag, y = prob_plot, fill = score)) +
  geom_text(data = source_totals,
            aes(x = flag, y = p_damaging,
                label = scales::percent(prop, accuracy = 0.1)),
            vjust = -2, size = 3.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  scale_x_discrete(expand = expansion(mult = c(1, 1))) + 
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks  = c(-5, 0, 5),
                       limits  = c(-5, 5),
                       labels  = c("benign", "unknown", "pathogenic"),
                       na.value = "grey80") +
  labs(x = NULL,
       y = "total\np(causal & damaging)") +
  theme_bw() +
  geom_text_repel(
    data = subset(variant_source_df, flag == "missing"),
    aes(x = flag, y = ymid, label = variant_lab, colour = score),
    nudge_x = -0.5,
    direction = "y",
    hjust = 1,
    size = 2,
    segment.size = 0.2,
    max.overlaps = 10,
    box.padding = 0.2
  ) +
  geom_text_repel(
    data = subset(variant_source_df, flag == "present"),
    aes(x = flag, y = ymid, label = variant_lab, colour = score),
    nudge_x = 0.5,
    direction = "y",
    hjust = 0,
    size = 2,
    segment.size = 0.2,
    max.overlaps = 10,
    box.padding = 0.2
  ) +
  scale_colour_gradientn(
    colours = c("navy", "blue", "darkred"),
    breaks  = c(-5, 0, 5),
    limits  = c(-5, 5),
    labels  = c("benign", "unknown", "pathogenic"),
    na.value = "grey80"
  ) +
  guides(color = "none") +
  geom_col(width = 0.3, colour = "black", position = position_stack())

# p_source_stack

# plot F: overall damaging‑causal distribution ----------------------------
# 1. posterior draws of relative share among causals
sim_causal_norm <- sim_causal / rowSums(sim_causal)   # n_sim × n_causal

# 2. build a matrix of P(G_i=1) for each causal variant i
#    known present → 1
#    known ref     → 0
#    missing       → use adj_occurrence_prob (p_i)
geno_prob <- risk_variants %>%
  filter(class == "causal") %>%
  transmute(variant_lab,
            g_prob = case_when(
              known_var_observed == 1 ~ 1,
              known_var_observed == 0 ~ 0,
              known_var_observed == -9 ~ adj_occurrence_prob
            )
  ) %>%
  pull(g_prob)
# geno_prob is a length‑n_causal vector in the same order as sim_causal_norm's columns

# 3. simulate G^(m) for each draw m
#    a matrix n_sim × n_causal where entry is Bernoulli(g_prob[i])
G_sim <- matrix(rbinom(n_sim * length(geno_prob), 1, rep(geno_prob, each = n_sim)),
                nrow = n_sim, ncol = length(geno_prob), byrow = FALSE)

# 4. for each draw, compute total probability that a causal is present
total_sim <- rowSums(sim_causal_norm * G_sim)

# 5. extract CI
ci_overall <- quantile(total_sim, c(0.025, 0.5, 0.975))

# 6. plot
p_overall <- ggplot(tibble(share = total_sim), aes(x = share)) +
  geom_histogram(bins = 50, fill = "red", alpha = 0.6, colour = "black") +
  geom_vline(xintercept = ci_overall[c(1, 3)], linetype = "dashed") +
  geom_vline(xintercept = ci_overall[2], colour = "purple") +
  annotate("text", 
           # x = ci_overall[2],
           x = 0.54,
           y = Inf,
           vjust = 1.1,
           label = sprintf("median = %.3f\n95%% CI = [%.3f, %.3f]",
                           ci_overall[2], ci_overall[1], ci_overall[3]),
           # hjust = 0.5, 
           size = 4) +
  labs(x = "total p(causal & damaging & present)", y = "count") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.8)))

# p_overall

# assemble panels ---------------------------------------------------------
p_quant <- (p_prior + p_dist + p_causal) /
  (p_source_stack + p_causal_dist) /
  p_overall +
  plot_layout(guides="collect",axis="collect",heights=c(2,1,1)) +
  plot_annotation(tag_levels="A", subtitle = paste0("Gene: ",  genetic_defect))

# p_quant
# ggsave("../images/plot_quant_uncert_ci.pdf", plot = p_quant, width = 9, height = 7)
ggsave(paste0("../images/plot_scenario_", scenario, "_quant_uncert_ci.pdf"), plot = p_quant, width = 9, height = 7)

# # tables ----
# library(knitr)
# library(kableExtra)

# first join in the CI columns (for present variants only)
variant_with_ci <- risk_variants %>%
  left_join(
    ci_present, by = "variant_lab"
  )

# 2) Build and round the table
variant_table_full <- variant_with_ci %>%
  mutate(
    lower  = as.numeric(lower),
    median = as.numeric(median),
    upper  = as.numeric(upper)
  ) %>%
  select(
    Variant         = variant_lab,
    Flag            = flag,
    Class           = class,
    Evidence_Score  = evidence_score,
    Occurrence_Prob = occurrence_prob,
    Adj_Occ_Prob    = adj_occurrence_prob,
    Alpha           = alpha,
    Beta            = beta,
    Lower           = lower,
    Median          = median,
    Upper           = upper,
    Posterior_Share = posterior_share,
    Prob_Causal     = prob_causal_damaging
  ) %>%
  # round to 3 significant digits, keep NAs as NA
  mutate(across(
    where(is.numeric),
    ~ if_else(is.na(.x), NA_real_, round(.x, 3))
  )) %>%
  arrange(desc(Prob_Causal))


# Add overall result to table
overall_row <- tibble(
  Variant         = "Total",
  Flag            = NA,
  Class           = NA,
  Evidence_Score  = NA,
  Occurrence_Prob = NA,
  Adj_Occ_Prob    = NA,
  Alpha           = NA,
  Beta            = NA,
  Lower           = round(ci_overall[1], 3),
  Median          = round(ci_overall[2], 3),
  Upper           = round(ci_overall[3], 3),
  Posterior_Share = NA,
  Prob_Causal     = round(ci_overall[2], 3),
)

# Bind to table
variant_table_full <- bind_rows(variant_table_full, overall_row)

# prep caption ----

# Identify top present and missing variants separately
top_present <- variant_table_full %>%
  filter(Flag == "present", !is.na(Prob_Causal), Variant != "Total") %>%
  arrange(desc(Prob_Causal)) %>%
  slice_head(n = 1)

top_missing <- variant_table_full %>%
  filter(Flag == "missing", !is.na(Prob_Causal), Variant != "Total") %>%
  arrange(desc(Prob_Causal)) %>%
  slice_head(n = 1)

# Extract the overall summary
overall_row <- variant_table_full %>%
  filter(Variant == "Total")

# Compose description fragments based on availability
desc_present <- if (nrow(top_present) > 0) {
  paste0("The most strongly supported observed variant was \\texttt{", 
         top_present$Variant, "} (posterior: ", top_present$Prob_Causal, "). ")
} else {
  "No observed variants were detected in this scenario. "
}

desc_missing <- if (nrow(top_missing) > 0) {
  paste0("The strongest unsequenced variant was \\texttt{", 
         top_missing$Variant, "} (posterior: ", top_missing$Prob_Causal, "). ")
} else {
  "No unsequenced variants were included in this scenario. "
}

# caption ----
caption_text <- paste0(
  "Result of clinical genetics diagnosis scenario ", scenario, " including metadata. ",
  desc_present,
  desc_missing,
  "The total probability of a causal diagnosis given the available evidence was ",
  overall_row$Prob_Causal, " (95\\% CI: ", overall_row$Lower, "--", overall_row$Upper, ").",
  "\\label{tab:table_scenario_", scenario, "_quant_uncert_ci}"
)

# caption = paste0("Result of clinical genetics diagnosis scenario ", scenario, ". The proband carried three observed variants, including the known pathogenic \\texttt{p.Ser237Ter} (true positive), and lacked coverage at three additional sites, including likely-pathogenic splice‑donor \\texttt{c.159+1G>A} (false negative). The damaging-only posterior probabilities for these two variants were 0.382 and 0.351, resulting in total probability (prob causal) of causal diagnosis given the existing evidence of 0.521 (95\\% CI: 0.248–0.787). \\label{tab:table_scenario_2_quant_uncert_ci}"),

# Generate the LaTeX tabuLower# Generate the LaTeX tabular code without the surrounding table environment
latex_tabular <- kable(
  variant_table_full,
  format      = "latex",
  linesep     = "",  # <- this removes all automatic row spacing
  booktabs    = TRUE,
  # caption = paste0("Result of clinical genetics diagnosis scenario ", scenario, ". The proband carried three observed variants, including the known pathogenic \\texttt{p.Ser237Ter} (true positive), and lacked coverage at three additional sites, including likely-pathogenic splice‑donor \\texttt{c.159+1G>A} (false negative). The damaging-only posterior probabilities for these two variants were 0.382 and 0.351, resulting in total probability (prob causal) of causal diagnosis given the existing evidence of 0.521 (95\\% CI: 0.248–0.787). \\label{tab:table_scenario_2_quant_uncert_ci}"),
  caption = caption_text,
  # escape      = FALSE,
  # table.envir = FALSE,                    # suppress \begin{table}...\end{table}
  col.names   = gsub("_", " ", names(variant_table_full))
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    position      = "center",
  ) %>%
  row_spec(6, extra_latex_after = "\\addlinespace") %>% 
  column_spec(1, width = "2.5cm") %>%  
  column_spec(2, width = "1.5cm") %>%  
  column_spec(3, width = "1.5cm") %>%  
  column_spec(4, width = "1.5cm") %>% 
  column_spec(5, width = "1.8cm") %>% 
  column_spec(6:10, width = "1.5cm") %>%
  column_spec(11:13, width = "1.5cm")

# Wrap in resizebox and full table environment
# cat(latex_tabular)

# write the LaTeX table to file
# out_file <- "../images/table_quant_uncert_ci.tex"
out_file <-  paste0("../images/table_scenario_", scenario, "_quant_uncert_ci.tex")
table_lines <- c(latex_tabular)
writeLines(table_lines, out_file)

# report ----
candidate_variants <- variant_table_full %>%
  filter(Class == "causal", Variant != "Total") %>%
  left_join(
    risk_variants %>%
      select(variant_lab, genename, Inheritance,
             HGVSc_VEP, HGVSp_VEP,
             gnomAD_genomes_AF, occurrence_prob,
             posterior_share, prob_causal_damaging, flag),
    by = c("Variant" = "variant_lab")
  ) %>%
  left_join(ci_present, by = c("Variant" = "variant_lab"))


# define parameters
param_names <- c(
  "Gene",
  "HGVSc",
  "HGVSp",
  "Inheritance",
  "Patient sex",
  "gnomAD frequency",
  "95% CI lower",
  "p(median)",
  "95% CI upper",
  "Posterior p(causal)",
  "Interpretation"
)

# build list keyed by flag ("present","missing")
vals_list <- lapply(seq_len(nrow(candidate_variants)), function(i) {
  row <- candidate_variants[i, ]
  c(
    Gene                   = row$genename,
    HGVSc                  = row$HGVSc_VEP,
    HGVSp                  = row$HGVSp_VEP,
    Inheritance            = row$Inheritance,
    `Patient sex`          = "Male",
    `gnomAD frequency`     = formatC(row$gnomAD_genomes_AF, format="e", digits=2),
    `95% CI lower`         = sprintf("%.3f", row$lower),
    `p(median)`            = sprintf("%.3f", row$median),
    `95% CI upper`         = sprintf("%.3f", row$upper),
    `Posterior p(causal)`  = sprintf("%.3f", row$prob_causal_damaging),
    Interpretation         = if (row$flag == "present") {
      "Reported causal; variant observed"
    } else {
      "Reported causal; variant not detected — consider follow‑up"
    }
  )
})
names(vals_list) <- candidate_variants$flag

# assemble wide table (columns = present, missing)
df_report_wide <- data.frame(
  Parameter = param_names,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
for(flag in unique(candidate_variants$flag)) {
  df_report_wide[[flag]] <- vals_list[[flag]][param_names]
}

# --- build caption including total summary ---
desc_present <- if(nrow(top_present)>0) {
  paste0("Reported causal: \\texttt{", top_present$Variant, "} (posterior ", top_present$Prob_Causal, "). ")
} else ""
desc_missing <- if(nrow(top_missing)>0) {
  paste0("Undetected causal: \\texttt{", top_missing$Variant, "} (posterior ", top_missing$Prob_Causal, "). ")
} else ""
desc_interp <- paste0(
  # desc_present, desc_missing,
  "Overall probability of correct causal diagnosis due to SNV/INDEL given the currently available evidence: ",
  overall_row$Prob_Causal,
  " (95\\% CI ", overall_row$Lower, "--", overall_row$Upper, "). "
)
# desc_interp

caption_report <- paste0(
  "Final variant report for clinical genetics scenario ", scenario, ". ",
  desc_present,
  desc_missing,
  "The total probability of a causal diagnosis given the available evidence was ",
  overall_row$Prob_Causal, " (95\\% CI: ", overall_row$Lower, "--", overall_row$Upper, ").",
  "\\label{tab:table_scenario_", scenario, "_report}"
)

n_cols <- ncol(df_report_wide)
# build a matching col.names vector:
#   first element is always "Parameter", then the flag names
col_names <- c("Parameter", setdiff(colnames(df_report_wide), "Parameter"))

# 1. Render the base table to character
latex_tab <- kable(
  df_report_wide,
  format    = "latex",
  linesep     = "",  # <- this removes all automatic row spacing
  booktabs  = TRUE,
  caption   = caption_report,
  col.names = c("Parameter", setdiff(names(df_report_wide), "Parameter"))
) %>%
  kable_styling(latex_options = c("hold_position","scale_down"),
                position = "center") %>%
  column_spec(1, width="4cm") %>%
  { 
    # widen the rest to match number of columns
    n_cols <- ncol(df_report_wide)
    column_spec(., 2:n_cols, width = paste0(12 / (n_cols-1), "cm")) 
  } %>%
  as.character()

# 2. Split into lines
lines <- strsplit(latex_tab, "\n", fixed = TRUE)[[1]]

# 3. Locate bottomrule
i_bot <- grep("^\\\\bottomrule", lines)

# 4. Build footer with “Summary” in col 1
n_cols <- ncol(df_report_wide)
footer <- paste0(
  "\\midrule\n",
  "\\textbf{Summary} & ",
  "\\multicolumn{", n_cols-1, "}{p{", paste0((n_cols-1)*5, "cm"), "}}{", 
  desc_interp, 
  "} \\\\"
)

new_lines <- append(lines, footer, after = i_bot - 1)

latex_tab_with_footer <- paste(new_lines, collapse = "\n")

out_file_report <-  paste0("../images/table_scenario_", scenario, "_final_report.tex")
table_lines_report <- c(latex_tab_with_footer)
writeLines(table_lines_report, out_file_report)
gc()
