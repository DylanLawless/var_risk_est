





# animate ----

# Start import data ----
header_line <- readLines("../data/nfkb1_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

# header_fields <- head(header_fields, 200)

library(ggplot2)
library(gganimate)
library(dplyr)
library(tibble)


# Create base data
base_df <- tibble(
  field = header_fields,
  value = seq_along(header_fields)
)

# Create cumulative data for animation
df <- lapply(seq_along(header_fields), function(i) {
  base_df %>%
    slice_head(n = i) %>%
    mutate(frame = i)
}) %>%
  bind_rows() %>%
  mutate(field = factor(field, levels = header_fields))

# Annotation label (top-left)
labels <- tibble(
  frame = seq_along(header_fields),
  label = paste("Evidence data:\n", header_fields)
)

# Plot
p <- ggplot(df, aes(x = value, y = value)) +
  geom_col(fill = "steelblue", color = "black") +
  geom_text(data = labels,
            aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.2, size = 5) +
  labs(x = "Database", y = "Cumaltive evidence sources") +
  coord_cartesian(ylim = c(0, length(header_fields) + 1)) +
  # theme_minimal(base_size = 12) +
  # theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.title = element_blank(),
    # panel.grid.minor = element_blank()
  # ) +
  transition_manual(frame)


anim0 <- animate(p, duration = 3, fps = 20, width = 700, height = 700, res = 120)
gif_filename <- file.path("./quant_anim0.gif")
anim_save(gif_filename, anim0)

# 1 ----
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)
library(scales)

# Setup for animation: interpolate x (occurrence_prob) from 0 to true value
n_frames <- 50

prior_df_anim <- prior_df %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(data = list(tibble(
    frame = 1:n_frames,
    x_interp = seq(0, occurrence_prob, length.out = n_frames),
    label_interp = sprintf("%.5f", seq(0, occurrence_prob, length.out = n_frames))
    # label_interp = formatC(seq(0, occurrence_prob, length.out = n_frames), format = "e", digits = 2)
  ))) %>%
  unnest(data)

p_prior_anim <- ggplot(prior_df_anim,
                       aes(x = x_interp, y = var_plot, fill = score, group = id)) +
  geom_point(shape = 21, colour = "black", size = 3, stroke = 0.3) +
  geom_text(aes(label = label_interp), hjust = 0, nudge_x = 0.0005, size = 3) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(group = label_wrap_gen(width = 9))) +
  # theme_minimal(base_size = 9) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks = c(-5, 0, 5),
                       labels = c("benign", "unknown", "pathogenic")) +
  labs(x = "prior occurrence probability", y = NULL,
       title = "Calculating probability occurrence priors") +
  scale_x_continuous(limits = c(0, 0.0035), breaks = pretty_breaks(n = 2)) +
  transition_time(frame) +
  ease_aes("cubic-in-out") 

# animate(p_prior_anim, fps = 20, duration = 3, width = 600, height = 400)
anim1 <- animate(p_prior_anim, duration = 3, fps = 20, width = 700, height = 700, res = 120)
gif_filename <- file.path("./quant_anim1.gif")
anim_save(gif_filename, anim1)


# 2 ----


library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)

# Parameters
n_points <- 1000
n_frames <- 30
bin_width <- 0.05
bins <- seq(0, 1, by = bin_width)

# Start: all zeros
start_values <- rep(0, n_points)

# End: clipped normal dist around 0.5
set.seed(42)
end_values <- pmin(pmax(rnorm(n_points, mean = 0.5, sd = 0.15), 0), 1)

# Interpolate values across frames
interp_df <- tibble(id = 1:n_points) %>%
  crossing(frame = 0:n_frames) %>%
  mutate(value = (frame / n_frames) * end_values)

# Bin values per frame
hist_df <- interp_df %>%
  mutate(bin = cut(value, breaks = bins, include.lowest = TRUE, labels = FALSE)) %>%
  count(frame, bin) %>%
  mutate(bin_mid = bins[bin] + bin_width / 2)


library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)
library(forcats)

# Parameters
n_frames <- 30
bin_width <- 0.05
bins <- seq(0, 1, by = bin_width)
bin_mids <- bins[-length(bins)] + bin_width / 2

# Get distinct variant info
variant_info <- share_df %>%
  distinct(variant_lab, var_plot, group, score)

# Build interpolated values per variant
interp_df <- share_df %>%
  group_by(variant_lab) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  crossing(frame = 0:n_frames) %>%
  mutate(value = (frame / n_frames) * share)

# Bin the interpolated values
hist_df <- interp_df %>%
  mutate(bin = cut(value, breaks = bins, include.lowest = TRUE, labels = FALSE)) %>%
  count(frame, variant_lab, bin, name = "count") %>%
  mutate(bin_mid = bin_mids[bin]) %>%
  left_join(variant_info, by = "variant_lab") %>%
  replace_na(list(count = 0))

# Combine group and var_plot into labelled facet strip
hist_df <- hist_df %>%
  mutate(var_lab = paste0(var_plot, "\n", group),
         var_lab = factor(var_lab, levels = unique(var_lab)))


medians <- interp_df %>%
  mutate(var_lab = paste0(var_plot, "\n", group)) %>%
  group_by(var_lab, frame) %>%
  summarise(med = median(value), .groups = "drop") %>%
  # mutate(label = paste0("median = ", formatC(med, format = "f", digits = 2)))
  mutate(label = paste0(formatC(med, format = "f", digits = 2)))


p_dist_anim <- ggplot(hist_df, aes(x = bin_mid, y = count, fill = score)) +
  geom_col(width = bin_width, colour = "black", alpha = 0.8) +
  facet_grid(var_lab ~ .,
             # scales = "free_y", 
             space = "free_y") +
  scale_fill_gradientn(
    colours = c("navy", "lightblue", "red"),
    breaks = c(-5, 0, 5),
    labels = c("benign", "unknown", "pathogenic")
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_text(data = medians,
            aes(x = med+0.2, y = Inf, label = label),
            hjust = 1, vjust = 1.2,
            size = 5, inherit.aes = FALSE) +
  labs(x = "posterior share", y = "count",
       title = "Calculating posterior share distribution (causal & damaging)") +
  # theme_minimal(base_size = 9) +
  theme(
    strip.text.y = element_text(size = 7),
    strip.placement = "outside",
    panel.spacing.y = unit(0.4, "lines")
  ) +
  transition_manual(frame)

# animate(p_dist_anim, fps = 20, duration = 3, width = 400, height = 400)
# anim2 <- animate(p_anim, nframes = 100, fps = 10, width = 700, height = 500, res = 120)
anim2 <- animate(p_dist_anim, duration = 3, fps = 20, width = 700, height = 700, res = 120)
gif_filename <- file.path("./quant_anim2.gif")
anim_save(gif_filename, anim2)

# ggsave(paste0("../images/plot_scenario_", scenario, "_quant_uncert_ci.pdf"), plot = p_quant, width = 9, height = 7)



# p_source_stack ----

library(ggplot2)
library(gganimate)
library(dplyr)
library(ggrepel)
library(scales)

# Base data
# variant_source_df <- tibble::tibble(
#   flag = c("missing", "missing", "missing", "present", "present", "present"),
#   score = c(-5, 0, 4.5, 0, 0, 5),
#   prob_causal_damaging = c(0, 0, 0.353, 0, 0, 0.381),
#   variant_lab = c("p.Thr567Ile", "p.Val236Ile", "c.159+1G>A", 
#                   "p.Arg231His", "p.Gly650Arg", "p.Ser237Ter"),
#   prob_plot = c(0.01, 0.01, 0.363, 0.01, 0.01, 0.391),
#   ymax = c(0.01, 0.02, 0.383, 0.01, 0.02, 0.411),
#   ymin = c(0, 0.01, 0.020, 0, 0.01, 0.020),
#   ymid = c(0.005, 0.015, 0.202, 0.005, 0.015, 0.215)
# )

# Dummy totals
source_totals <- variant_source_df %>%
  group_by(flag) %>%
  summarise(p_damaging = sum(prob_plot), .groups = "drop") %>%
  mutate(prop = p_damaging / sum(p_damaging))

# Create interpolated frames
n_frames <- 30
interp_df <- variant_source_df %>%
  slice(rep(1:n(), each = n_frames)) %>%
  group_by(variant_lab) %>%
  mutate(frame = 1:n_frames,
         prob_plot_interp = (frame / n_frames) * prob_plot,
         ymid_interp = (frame / n_frames) * ymid)

library(scales)

# Create interpolated source totals
n_frames <- 30
animated_totals <- source_totals %>%
  slice(rep(1:n(), each = n_frames)) %>%
  group_by(flag) %>%
  mutate(
    frame = 1:n_frames,
    p_damaging_interp = (frame / n_frames) * p_damaging,
    prop_interp = (frame / n_frames) * prop,
    label = percent(prop_interp, accuracy = 0.1)
  )

# Animation
p_source_stack_anim <- ggplot(interp_df, aes(x = flag, y = prob_plot_interp, fill = score)) +
  geom_col(width = 0.3, colour = "black", position = "stack") +
  # geom_text(data = source_totals,
  #           aes(x = flag, y = p_damaging,
  #               label = percent(prop, accuracy = 0.1)),
  #           vjust = -2, size = 3.5, inherit.aes = FALSE) +
  geom_text(data = animated_totals,
            aes(x = flag, y = p_damaging_interp, label = label),
            vjust = -2, size = 3.5, inherit.aes = FALSE) +
  geom_text_repel(
    data = subset(interp_df, flag == "missing"),
    aes(x = flag, y = ymid_interp, label = variant_lab, colour = score),
    nudge_x = -0.5, direction = "y", hjust = 1, size = 2,
    segment.size = 0.2, max.overlaps = 10, box.padding = 0.2
  ) +
  geom_text_repel(
    data = subset(interp_df, flag == "present"),
    aes(x = flag, y = ymid_interp, label = variant_lab, colour = score),
    nudge_x = 0.5, direction = "y", hjust = 0, size = 2,
    segment.size = 0.2, max.overlaps = 10, box.padding = 0.2
  ) +
  scale_fill_gradientn(
    colours = c("navy", "lightblue", "red"),
    breaks = c(-5, 0, 5), limits = c(-5, 5),
    labels = c("benign", "unknown", "pathogenic"), na.value = "grey80"
  ) +
  scale_colour_gradientn(
    colours = c("navy", "blue", "darkred"),
    breaks = c(-5, 0, 5), limits = c(-5, 5),
    labels = c("benign", "unknown", "pathogenic"), na.value = "grey80"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  scale_x_discrete(expand = expansion(mult = c(1, 1))) +
  labs(x = NULL, y = "total p(causal & damaging)",
       title = "Calculating posterior total share of causal & damaging") +
  theme_bw() +
  guides(color = "none") +
  transition_manual(frame)

# Animate
# animate(p_anim, fps = 20, width = 800, height = 500, res = 120)
anim3 <- animate(p_source_stack_anim, duration = 3, fps = 20, width = 700, height = 700, res = 120)
gif_filename <- file.path("./quant_anim3.gif")
anim_save(gif_filename, anim3)





# p_anim_dist ----

library(ggplot2)
library(gganimate)
library(dplyr)
library(tibble)

# 1. use the calculated probabilities
set.seed(1)
# total_sim <- rowSums(matrix(runif(1000), ncol = 1))  # replace with your real `total_sim`

# split the data into two groups
n <- length(total_sim)
sim_df <- tibble(
  id = 1:n,
  original = total_sim,
  group = sample(c("to_0", "to_1"), size = n, replace = TRUE)
)

# interpolate over frames
n_frames <- 30
anim_df <- sim_df %>%
  slice(rep(1:n(), each = n_frames)) %>%
  group_by(id) %>%
  mutate(
    frame = 1:n_frames,
    value = case_when(
      group == "to_0" ~ original * (1 - frame / n_frames),
      group == "to_1" ~ original * (1 - frame / n_frames) + (frame / n_frames)
    )
  )


# Define colour ramps
ramp_to_0 <- colorRampPalette(c("pink", "blue"))
ramp_to_1 <- colorRampPalette(c("pink", "red"))

# Add frame-specific colour per row
anim_df <- anim_df %>%
  mutate(
    colour = case_when(
      group == "to_0" ~ ramp_to_0(n_frames)[frame],
      group == "to_1" ~ ramp_to_1(n_frames)[frame]
    )
  )

median_labels <- anim_df %>%
  group_by(group, frame) %>%
  summarise(med = median(value), .groups = "drop") %>%
  mutate(
    label = case_when(
      group == "to_0" ~ sprintf("Resolved missing\nCounterfactual: %.2f\n", med),
      group == "to_1" ~ sprintf("Present variant\nCausal: %.2f\n", med)
    )
  )

median_labels <- anim_df %>%
  group_by(group, frame) %>%
  summarise(
    q025 = quantile(value, 0.025),
    med  = median(value),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      group == "to_0" ~ sprintf(
        "Resolved missing\n95%% CI: [%.2f, %.2f]\nCounterfactual: %.2f", q025, q975, med
      ),
      group == "to_1" ~ sprintf(
        "Present variant\n95%% CI: [%.2f, %.2f] \nCausal: %.2f", q025, q975, med
      )
    )
  )


p_anim_dist <- ggplot(anim_df, aes(x = value)) +
  # geom_histogram(bins = 50, fill = "red", alpha = 0.6, colour = "black") +
  geom_histogram(aes(fill = colour), bins = 50, alpha = 0.6, colour = "black") +
  scale_fill_identity() +
  geom_text(data = median_labels,
            aes(x = med, y = (med*500)+500, label = label),
            inherit.aes = FALSE,
            vjust = -0.5, hjust = 0.5, size = 4, colour = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1500), 
                     # expand = expansion(mult = c(0, 0.1))
                     ) +
  labs(
    title = "Resolving total p(causal & damaging & present)",
    x = "p(causal & damaging & present)",
    y = "count"
  ) +
  # theme_minimal(base_size = 12) +
  transition_manual(frame)

# animate(p_anim_dist, fps = 20, width = 700, height = 400, res = 120)

anim4 <- animate(p_anim_dist, duration = 3, fps = 20, width = 700, height = 700, res = 120)
gif_filename <- file.path("./quant_anim4.gif")
anim_save(gif_filename, anim4)



