
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
sim_mat <- sapply(seq_len(nrow(risk_variants)),
                  function(i) rbeta(n_sim,
                                    risk_variants$alpha[i],
                                    risk_variants$beta[i]))
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
prior_df <- risk_variants %>%
  mutate(var_plot = reorder_within(variant_lab, occurrence_prob, group))

p_prior <- ggplot(prior_df,
                  aes(x = occurrence_prob, y = var_plot, fill = score)) +
  geom_point(shape = 21, colour = "black", size = 3, stroke = 0.3) +
  
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks = c(-5, 0, 5),
                       labels = c("benign", "unknown", "pathogenic")) +
  scale_y_reordered() +
  labs(x = "prior occurrence probability", y = NULL) +
  scale_x_continuous(limits = c(0, 0.003),
                     breaks = pretty_breaks(n = 2))

# — plot B: full posterior distributions ——————————————————————————
share_df <- as_tibble(share_all) %>%
  mutate(draw=row_number()) %>%
  pivot_longer(-draw, names_to="variant_lab", values_to="share") %>%
  left_join(risk_variants %>% select(variant_lab, group, score),
            by="variant_lab") %>%
  mutate(var_plot=reorder_within(variant_lab, share, group))

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
p_causal <- ggplot(
  risk_variants %>% mutate(var_plot=reorder_within(
    variant_lab, prob_causal_damaging, group)),
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
  select(flag, score, prob_causal_damaging) %>%
  mutate(
    # add epsilon so zero‑prob variants still get a sliver
    prob_plot = prob_causal_damaging + eps
  )

source_totals <- variant_source_df %>%
  group_by(flag) %>%
  summarise(p_damaging = sum(prob_causal_damaging), .groups = "drop") %>%
  mutate(prop = p_damaging / sum(p_damaging))

p_source_stack <- ggplot(variant_source_df,
                         aes(x = flag, y = prob_plot, fill = score)) +
  geom_col(width = 0.6, colour = "black", position = position_stack()) +
  geom_text(data = source_totals,
            aes(x = flag, y = p_damaging,
                label = scales::percent(prop, accuracy = 0.1)),
            vjust = -2, size = 3.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks  = c(-5, 0, 5),
                       limits  = c(-5, 5),
                       labels  = c("benign", "unknown", "pathogenic"),
                       na.value = "grey80") +
  labs(x = NULL,
       y = "total\np(causal & damaging)") +
  theme_bw()

p_source_stack


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
           x = 0.5,
           y = Inf,
           vjust = 1.2,
           label = sprintf("median = %.3f\n95%% CI = [%.3f, %.3f]",
                           ci_overall[2], ci_overall[1], ci_overall[3]),
           hjust = 0.5, size = 4) +
  labs(x = "total p(causal & damaging & present)", y = "count") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.6)))

p_overall

# assemble panels ---------------------------------------------------------
p_quant <- (p_prior + p_dist + p_causal) /
  (p_source_stack + p_causal_dist) /
  p_overall +
  plot_layout(guides="collect",axis="collect",heights=c(2,1,1)) +
  plot_annotation(tag_levels="A")

p_quant
# ggsave("../images/plot_quant_uncert_ci.pdf", plot = p_quant, width = 9, height = 7)
ggsave(paste0("../images/plot_", scenario, "_quant_uncert_ci.pdf"), plot = p_quant, width = 9, height = 7)

# tables ----
library(knitr)
library(kableExtra)

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

class(variant_table_full$Lower)
variant_table_full$Lower

# Generate the LaTeX tabuLower# Generate the LaTeX tabular code without the surrounding table environment
latex_tabular <- kable(
  variant_table_full,
  format      = "latex",
  booktabs    = TRUE,
  # escape      = FALSE,
  # table.envir = FALSE,                    # suppress \begin{table}...\end{table}
  col.names   = gsub("_", " ", names(variant_table_full))
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    position      = "center"
  ) %>%
  column_spec(1, width = "2.5cm") %>%  
  column_spec(2, width = "1.8cm") %>%  
  column_spec(3, width = "1.8cm") %>%  
  column_spec(4, width = "1.5cm") %>% 
  column_spec(5, width = "1.8cm") %>% 
  column_spec(6:10, width = "1.5cm") %>%
  column_spec(11:13, width = "1.8cm")

# Wrap in resizebox and full table environment
cat(latex_tabular)

# write the LaTeX table to file
# out_file <- "../images/table_quant_uncert_ci.tex"
out_file <-  paste0("../images/table_", scenario, "_quant_uncert_ci.tex")

table_lines <- c(latex_tabular)

writeLines(table_lines, out_file)
