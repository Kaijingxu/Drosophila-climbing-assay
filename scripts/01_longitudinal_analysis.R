library(tidyverse)
library(lme4)
library(emmeans)
library(readxl)
library(scales)

# Import and clean raw longitudinal climbing data 
df <- readxl::read_excel("data/longitudinal.xlsx") %>%
  tidyr::fill(`Fly No.`, Genotype, Age, Round, Tray, .direction = "down") %>% 
  rename(fly_id_raw = `Fly No.`,
           genotype = Genotype,
          trial_num = Triplicate,
           age_day  = Age,
         round_id   = Round,
         tray_local = Tray) %>%
  mutate(fly_id    = factor(fly_id_raw),
         round_id  = factor(round_id),
         tray_local= factor(tray_local),
         trial_num = factor(trial_num),
       age_day_num = readr::parse_number(age_day),
       age_fac     = factor(age_day_num,
                            levels = c(1,7,14,21,28),
                            labels = paste0("Day ", c(1,7,14,21,28)),ordered = TRUE),
       genotype = factor(genotype, levels = c("+/+", "s64w/+", "s64w/s64w")))


# Convert interval scores to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("Interval"),
               names_to = "interval",
               names_prefix = "Interval ",
               values_to = "success",
               names_transform = list(interval = readr::parse_number)) %>%
  mutate(interval  = factor(interval, levels = c(1, 2, 3, 4), ordered = TRUE),
         interval_num = as.integer(interval),
         success      = as.integer(success)) %>%
  filter(!is.na(success), 
         success %in% c(0L, 1L), 
         !is.na(genotype),
         !is.na(age_fac))

# A graph: longitudinal baseline GLMM-predicted success
#GLMM model (formulation)
m_water <- glmer(success ~ genotype * age_fac + genotype * poly(interval_num, 2) + 
                   trial_num +
                   (1 | fly_id) +
                   (1 | age_fac:round_id) +
                   (1 | age_fac:round_id:tray_local),
                 data = df_long,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5)))
summary(m_water)

# Pearson overdispersion check
overdisp <- sum(resid(m_water, type = "pearson")^2) / df.residual(m_water)
cat("Overdispersion ratio =", overdisp, "\n")

# GLMM-predicted response probabilities by genotype, age and interval
emm_age <- emmeans(m_water,
                   ~ genotype | age_fac * interval_num,
                   at = list(interval_num = 1:4),
                   weights = "proportional")
  
emm_age_resp <- regrid(emm_age, transform = "response")

plot_data_age <- summary(emm_age_resp, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(interval_num = as.integer(interval_num),
         genotype     = fct_relevel(genotype, "+/+","s64w/+","s64w/s64w"))

# Contrast: homozygous mutant versus pooled wild type and heterozygote
genotype_levels <- levels(emm_age_resp)$genotype
mut_vs_pooled <- c("+/+" = -0.5,
                   "s64w/+" = -0.5,
                   "s64w/s64w" = 1)

mut_vs_pooled <- setNames(mut_vs_pooled[genotype_levels], genotype_levels)

cmp_age <- contrast(emm_age_resp,
                    method = list("s64w/s64w − pooled (+/+ and s64w/+)" = mut_vs_pooled),
                    by     = c("age_fac", "interval_num"),
                    adjust = "none") 

out_rd <- summary(cmp_age, infer = TRUE) %>%
  as.data.frame() %>%
  group_by(age_fac) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm")) %>%
  ungroup() %>%
  mutate(RD_pct  = estimate * 100,
         LCL_pct = asymp.LCL * 100,
         UCL_pct = asymp.UCL * 100,
         label   = if_else(p_holm < 1e-3, sprintf("p=%.1e", p_holm),
                          paste0("p=", formatC(p_holm, format = "f", digits = 3)))) %>%
  select(age_fac, interval_num, contrast, estimate, SE, asymp.LCL, asymp.UCL,
         RD_pct, LCL_pct, UCL_pct, z.ratio, p.value, p_holm, label)

# Label positions for p-values
label_positions <- plot_data_age %>%
  filter(genotype == "s64w/s64w") %>%
  transmute(age_fac, 
            interval_num, 
            prob, 
            asymp.LCL, 
            asymp.UCL,
            headroom = 1 - asymp.UCL,
            y = case_when(headroom >= 0.03 ~ pmin(0.98, asymp.UCL + 0.01),
                          prob >= 0.12     ~ pmax(0.02, prob - 0.10),
                          TRUE             ~ pmin(0.98, prob + 0.01))) %>%
  left_join(out_rd %>% transmute(age_fac, interval_num, label),
            by = c("age_fac", "interval_num")) %>%
  filter(!is.na(label))

# Fly numbers per genotype and day for legend labels 
n_by_day <- df_long %>%
  distinct(age_fac, genotype, fly_id) %>%
  count(age_fac, genotype, name = "n_day")

n_summary <- n_by_day %>%
  group_by(genotype) %>%
  summarise(n_min = min(n_day), n_max = max(n_day), .groups = "drop") %>%
  mutate(lbl = if_else(n_min == n_max,
                       paste0(as.character(genotype), 
                              " (n=", n_min, "/day)"),
                       paste0(as.character(genotype), 
                              " (n=", n_min, "–", n_max, "/day)")))

lab_vec <- setNames(n_summary$lbl, as.character(n_summary$genotype))


# order age facets
age_levels <- paste0("Day ", c(1,7,14,21,28))

plot_data_age <- plot_data_age %>%
  mutate(age_fac = factor(age_fac, levels = age_levels, ordered = TRUE))

out_rd <- out_rd %>%
  mutate(age_fac = factor(age_fac, levels = age_levels, ordered = TRUE))

label_positions <- label_positions %>%
  mutate(age_fac = factor(age_fac, levels = age_levels, ordered = TRUE))

# plot panel A
pd <- position_dodge(width = 0.20)

p_figA <- ggplot(plot_data_age,
                 aes(x = interval_num, y = prob, color = genotype, shape = genotype, group = genotype)) +
  geom_line(linewidth = 1, position = pd) +
  geom_point(size = 2.6, position = pd) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.10, position = pd, show.legend = FALSE) +
  facet_wrap(~ age_fac, nrow = 1) +
  scale_x_continuous(breaks = 1:4, labels = paste("Interval", 1:4)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_color_manual(name   = "Genotype",
                     values = c("+/+" = "blue", "s64w/+" = "green", "s64w/s64w" = "red"),
                     labels = lab_vec[levels(plot_data_age$genotype)]) +
  scale_shape_manual(name   = "Genotype",
                     values = c("+/+" = 16, "s64w/+" = 17, "s64w/s64w" = 15)) +
  guides(shape = "none",
         color = guide_legend(override.aes = list(shape = c(16, 17, 15), linewidth = 1, size = 2.6))) +
  labs(title = "A: Baseline GLMM-predicted success by test day and interval",
       x = "Interval", y = "Predicted probability of success", color = "Genotype") +
  theme_bw() +
  theme(strip.background    = element_rect(fill = "white", colour = "black"),
        panel.grid.minor    = element_blank(),
        plot.title.position = "plot",
        plot.title          = element_text(hjust = 0, face = "bold", size = 16),
        legend.title        = element_text(face = "bold")) +
  geom_text(data = label_positions, aes(x = interval_num, y = y, label = label),
            inherit.aes = FALSE, size = 3, vjust = 0)

p_figA

# Panel B: genotype deficit pooled across intervals
# Δ = 0.5 * P(+/+) + 0.5 * P(s64w/+) - P(s64w/s64w)
emm_age_avg <- emmeans(m_water,
                       specs = ~ genotype | age_fac,
                       at = list(interval_num = 1:4),
                       weights = "proportional")

emm_age_avg_resp <- regrid(emm_age_avg, transform = "response")

genotype_levels <- levels(emm_age_avg_resp)$genotype

deficit_contrast <- c("+/+" = 0.5,
                      "s64w/+" = 0.5,
                      "s64w/s64w" = -1)

deficit_contrast <- setNames(deficit_contrast[genotype_levels],
                             genotype_levels)

cmp_deficit <- contrast(emm_age_avg_resp,
                        method = list("pooled(+/+ and s64w/+) - s64w/s64w" = deficit_contrast),
                        by = "age_fac",
                        adjust = "none")

plot_data_deficit <- summary(cmp_deficit, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(age_fac = factor(age_fac,
                          levels = paste0("Day ", c(1, 7, 14, 21, 28)),
                          ordered = TRUE),
         deficit_pp = estimate * 100,
         LCL_pp = asymp.LCL * 100,
         UCL_pp = asymp.UCL * 100) %>%
  arrange(age_fac) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm"),
         p_label = if_else(p_holm < 1e-3,
                           sprintf("p=%.1e", p_holm),
                           paste0("p=", formatC(p_holm, digits = 3, format = "f"))),
         y_label = if_else(deficit_pp >= 0, UCL_pp + 2, LCL_pp - 2))

day21_x <- which(levels(plot_data_deficit$age_fac) == "Day 21")

#plot panel B
p_figB <- ggplot(plot_data_deficit, aes(x = age_fac, y = deficit_pp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect",
           xmin = day21_x - 0.25,
           xmax = day21_x + 0.25,
           ymin = -Inf,
           ymax = Inf,
           alpha = 0.08) +
  geom_pointrange(aes(ymin = LCL_pp, ymax = UCL_pp)) +
  geom_text(aes(y = y_label, label = p_label), vjust = 0, size = 4) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = " pp"),
                     expand = expansion(mult = c(0.02, 0.14))) +
  coord_cartesian(clip = "off") +
  labs(title = "B: Δ pooled across intervals — Δ = 0.5·P(WT) + 0.5·P(Het) − P(mut)",
       x = "Test day",
       y = expression(Delta~"success probability (pp)")) +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(face = "bold"))

p_figB

# Panel C: floor occupancy in homozygous mutants
# Floor occupancy = percentage of mutant flies with 0/3 successes
mutant_floor_by_fly <- df_long %>%
  filter(genotype == "s64w/s64w") %>%
  group_by(age_fac, interval_num, fly_id) %>%
  summarise(all_zero = as.integer(sum(success) == 0),
            n_trials = n(),
            .groups = "drop") %>%
  filter(n_trials >= 1)

wilson_ci <- function(x, n, conf = 0.95) 
  {if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p_hat <- x / n
  denominator <- 1 + z^2 / n
  centre <- p_hat + z^2 / (2 * n)
  adjustment <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n)
  c(lower = (centre - adjustment) / denominator,
    upper = (centre + adjustment) / denominator)}

floor_summary <- mutant_floor_by_fly %>%
  group_by(age_fac, interval_num) %>%
  summarise(n_fly = n(),
            n_floor = sum(all_zero),
            .groups = "drop") %>%
  group_by(age_fac) %>%
  summarise(n_fly_total = sum(n_fly),
            n_floor_total = sum(n_floor),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(floor_prop = n_floor_total / n_fly_total,
         ci = list(wilson_ci(n_floor_total, n_fly_total)),
         LCL = ci[1],
         UCL = ci[2]) %>%
  ungroup() %>%
  mutate(age_fac = factor(age_fac,
                          levels = paste0("Day ", c(1, 7, 14, 21, 28)),
                          ordered = TRUE),
         floor_pct = 100 * floor_prop,
         LCL_pct = 100 * LCL,
         UCL_pct = 100 * UCL,
         is_selected = age_fac == "Day 21")

day21_x <- which(levels(floor_summary$age_fac) == "Day 21")
label_padding <- 3
y_top <- min(100, max(floor_summary$UCL_pct, na.rm = TRUE) + label_padding)

#plot panel C
p_figC <- ggplot(floor_summary, aes(x = age_fac, y = floor_pct)) +
  annotate("rect",  xmin = day21_x - 0.45,
           xmax = day21_x + 0.45,
           ymin = -Inf,
           ymax = Inf,
           alpha = 0.08) +
  geom_col(aes(fill = is_selected), width = 0.65) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "grey40"),
                    guide = "none") +
  geom_errorbar(aes(ymin = LCL_pct, ymax = UCL_pct), width = 0.15) +
  geom_text(aes(y = UCL_pct + label_padding, label = sprintf("%.0f%%", floor_pct)),
            vjust = 0,
            size = 4) +
  coord_cartesian(ylim = c(0, y_top), clip = "off") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "%"),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "C: Floor occupancy in homozygous mutants",
       x = "Test day",
       y = "Floor occupancy across intervals (%)") +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

p_figC