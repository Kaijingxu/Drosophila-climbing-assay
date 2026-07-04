library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(DHARMa)
library(scales)
library(patchwork)
library(binom)
library(performance)
library(writexl)

# Salbutamol worked application: data import and model fitting
df <- readxl::read_excel("data/salbutamol_application.xlsx") %>%
  fill(`Fly No.`, Genotype, Conditions, Round, Tray, .direction = "down") %>%
  rename(fly_id    = `Fly No.`,
         genotype  = Genotype,
         trial_num = Triplicate,
         Condition = Conditions,
         round_id  = Round,
         tray_local = Tray) %>%
  pivot_longer(cols = starts_with("Interval"),
               names_to = "interval",
               names_prefix = "Interval ",
               values_to = "success") %>%
  filter(!is.na(success)) %>%
  mutate(fly_id       = factor(fly_id),
         genotype     = fct_relevel(factor(genotype), "+/+", "s64w/+", "s64w/s64w"),
         Condition    = fct_relevel(factor(Condition), "DMSO", "Sal"),
         trial_num    = factor(trial_num),
         round_id     = factor(round_id),
         tray_local   = factor(tray_local),
         interval     = factor(interval, levels = c("1", "2", "3", "4"), ordered = TRUE),
         interval_num = as.numeric(as.character(interval)),
         success      = as.integer(success))

# GLMM model (formulation)
m <- glmer(success ~ genotype * Condition * poly(interval_num, 2) + trial_num +
                 (1 | fly_id) +
                 (1 | round_id) +
                 (1 | round_id:tray_local),
               data = df,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))

summary(m)

# Pearson overdispersion check
overdisp <- sum(resid(m, type = "pearson")^2) / df.residual(m)
cat("Overdispersion ratio =", overdisp, "\n")


# GLMM-predicted probabilities by genotype, condition and interval
emm <- emmeans(m, ~ Condition | genotype * interval_num, 
               at = list(interval_num = c(1,2,3,4)),
               weights = "proportional")

plot_data <- summary(emm, type = "response", infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(interval_num = as.numeric(as.character(interval_num)))

# Response-scale Sal - DMSO contrasts by genotype and interval
emm_resp <- regrid(emm, transform = "response")

contrast_coef <- c("DMSO" = -1, "Sal" = 1)

cmp_df <- contrast(emm_resp,
                   method = list("Sal - DMSO" = contrast_coef),
                   by = c("genotype", "interval_num"),
                   adjust = "none") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(interval_num = as.numeric(as.character(interval_num)),
         RD_pct = estimate * 100,
         LCL_pct = asymp.LCL * 100,
         UCL_pct = asymp.UCL * 100) %>%
  group_by(genotype) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm"),
         label = case_when(
           p_holm < 1e-6  ~ sprintf("p = %.1e", p_holm),
           p_holm > 0.999 ~ "p > 0.999",
           TRUE           ~ paste0("p = ", formatC(p_holm, digits = 3, format = "f")))) %>%
  ungroup()

# P-value label positions
pos_df <- plot_data %>%
  group_by(genotype, interval_num) %>%
  summarise(y = pmin(1.05, max(asymp.UCL, na.rm = TRUE) + 0.04),
            .groups = "drop") %>%
  left_join(cmp_df %>% select(genotype, interval_num, label),
            by = c("genotype", "interval_num"))

# Optional: fly numbers for figure captions or supplementary information
fly_n <- df %>%
  group_by(genotype, Condition) %>%
  summarise(n_fly = n_distinct(fly_id), .groups = "drop") %>%
  pivot_wider(names_from = Condition,
              values_from = n_fly,
              names_prefix = "n_")
cap <- fly_n %>%
  transmute(txt = sprintf("%s: flies DMSO=%d, Sal=%d",
                          genotype, n_DMSO, n_Sal)) %>%
  pull(txt) %>%
  paste(collapse = " | ")

# Main figure: predicted probabilities
p_prob_panelA <- ggplot(plot_data,
                 aes(x = interval_num, y = prob,
                     color = Condition,
                     linetype = Condition,
                     shape = Condition)) +
  geom_point(position = position_dodge(width = 0.2), size = 2.5) +
  geom_line(position = position_dodge(width = 0.2), linewidth = 1) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1,
                position = position_dodge(width = 0.2)) +
  geom_text(data = pos_df,
            aes(x = interval_num, y = y, label = label),
            inherit.aes = FALSE,
            size = 3) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_x_continuous(breaks = 1:4,
                     labels = paste("Interval", 1:4)) +
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_color_manual(values = c("DMSO" = "grey30", "Sal" = "orange2")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = "Interval",
       y = "Predicted probability of success",
       color = "Condition",
       linetype = "Condition",
       shape = "Condition",
       title = "A: GLMM-predicted climbing success in the Day 21 vehicle-salbutamol application") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(face = "bold"),
        plot.margin = margin(t = 16, r = 8, b = 8, l = 8))

p_prob_panelA

# Interval-pooled Sal - DMSO effect by genotype
# The emmeans interaction note is expected because Condition interacts with interval.
emm_overall <- emmeans(m, ~ Condition | genotype,
                       at = list(interval_num = 1:4),
                       weights = "proportional")

emm_resp_overall <- regrid(emm_overall, transform = "response")

# Δprob=p^​(Sal)−p^​(DMSO)
delta_overall <- contrast(emm_resp_overall,
                          method = "revpairwise",
                          by = "genotype",
                          adjust = "none") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  rename(delta_prob = estimate,
         SE_delta   = SE,
         p_raw      = p.value) %>%
  mutate(genotype = factor(genotype, levels = c("+/+", "s64w/+", "s64w/s64w")),
         p_holm = p.adjust(p_raw, method = "holm"),
         p_label = case_when(p_holm < 1e-6  ~ sprintf("p = %.1e", p_holm),
                             p_holm > 0.999 ~ "p > 0.999",
                             TRUE           ~ paste0("p = ", formatC(p_holm, digits = 3, format = "f"))))

# plot
p_delta_panelB <- ggplot(delta_overall,
                         aes(x = genotype, y = delta_prob, fill = genotype)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_col(width = 0.55, colour = "black") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1, linewidth = 0.6) +
  geom_text(aes(y = asymp.UCL + 0.04, label = p_label),
            vjust = 0, size = 5) +
  scale_y_continuous(labels = function(x) paste0(x * 100, " pp"),
                     breaks = seq(-0.2, 1, by = 0.2),
                     expand = expansion(mult = c(0.02, 0.12))) +
  coord_cartesian(ylim = c(-0.2, 1), clip = "off") +
  scale_fill_manual(values = c("+/+" = "blue",
                                "s64w/+" = "green",
                                "s64w/s64w" = "red")) +
  labs(x = "Genotype",
       y = "Δ success probability (Sal - DMSO, pp)",
       fill = "Genotype",
       title = "B: Interval-pooled salbutamol effect by genotype") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(face = "bold"))

p_delta_panelB

# Optional: between-genotype contrast of treatment effects (ΔΔp)
delta_overall_emm <- contrast(emm_resp_overall,
                              method = "revpairwise",
                              by = "genotype",
                              adjust = "none")

genotype_levels <- levels(delta_overall_emm)$genotype

ddp_weights <- c("+/+" = -0.5,
                 "s64w/+" = -0.5,
                 "s64w/s64w" = 1)

ddp_weights <- setNames(ddp_weights[genotype_levels], genotype_levels)

ddp_out <- contrast(delta_overall_emm,
                    method = list("mutant - pooled WT/Het" = ddp_weights),
                    by = "contrast") %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# Optional: subtitle for p_delta_panelB
p_txt <- ifelse(ddp_out$p.value < 1e-6,
                sprintf("p = %.1e", ddp_out$p.value),
                paste0("p = ", formatC(ddp_out$p.value, digits = 3, format = "f")))

subtitle_line <- sprintf("ΔΔp = %.1f%% (95%% CI %.1f–%.1f%%), %s",
                         100 * ddp_out$estimate,
                         100 * ddp_out$asymp.LCL,
                         100 * ddp_out$asymp.UCL,
                         p_txt)
# p_delta_panelB <- p_delta_panelB + labs(subtitle = subtitle_line)


# Round-level sensitivity analyses
agg_int <- df %>%
  group_by(genotype, round_id, Condition, interval_num) %>%
  summarise(success = sum(success),
            total = n(),
            .groups = "drop") %>%
  mutate(prop = (success + 0.5) / (total + 1),
         var  = prop * (1 - prop) / (total + 1),
         w    = 1 / pmax(var, 1e-6))

B <- 9999
ALT <- "greater"

perm_vec <- function(x_Sal, x_DMSO, B = 9999,
                     alternative = c("two.sided", "less", "greater"),
                     seed = 123) {
  alternative <- match.arg(alternative)
  
  if (length(x_Sal) == 0L || length(x_DMSO) == 0L) {
    return(NA_real_)
  }
  
  set.seed(seed)
  
  x <- c(x_Sal, x_DMSO)
  n1 <- length(x_Sal)
  obs <- mean(x_Sal) - mean(x_DMSO)
  
  reps <- replicate(B, {
    idx <- sample(seq_along(x), n1, replace = FALSE)
    mean(x[idx]) - mean(x[-idx])
  })
  
  if (alternative == "greater") {
    (sum(reps >= obs) + 1) / (B + 1)
  } else if (alternative == "less") {
    (sum(reps <= obs) + 1) / (B + 1)
  } else {
    (sum(abs(reps) >= abs(obs)) + 1) / (B + 1)
  }
}

perm_panel_fast <- function(y, g, B = 9999, alt = c("greater", "less")) {
  alt <- match.arg(alt)
  n1 <- sum(g == "Sal")
  n <- length(y)
  
  if (n1 == 0 || n1 == n) {
    return(c(stat = NA_real_, p = NA_real_))
  }
  
  stat_obs <- mean(y[g == "Sal"]) - mean(y[g == "DMSO"])
  idx_mat <- replicate(B, sample.int(n, n1))
  z1 <- colMeans(matrix(y[idx_mat], nrow = n1))
  z0 <- (sum(y) - colSums(matrix(y[idx_mat], nrow = n1))) / (n - n1)
  reps <- z1 - z0
  
  p <- if (alt == "greater") {
    (sum(reps >= stat_obs) + 1) / (B + 1)
  } else {
    (sum(reps <= stat_obs) + 1) / (B + 1)
  }
  
  c(stat = stat_obs, p = p)
}

# Optional: per-interval round-level permutation sensitivity
sal_interval_table <- agg_int %>%
  select(genotype, interval_num, round_id, Condition, prop) %>%
  pivot_wider(names_from = Condition, values_from = prop) %>%
  group_by(genotype, interval_num) %>%
  summarise(n_round_Sal  = sum(!is.na(Sal)),
            n_round_DMSO = sum(!is.na(DMSO)),
            mean_Sal     = mean(Sal, na.rm = TRUE),
            mean_DMSO    = mean(DMSO, na.rm = TRUE),
            delta_prob   = mean_Sal - mean_DMSO,
            vec_Sal      = list(Sal[!is.na(Sal)]),
            vec_DMSO     = list(DMSO[!is.na(DMSO)]),
            .groups = "drop") %>%
  mutate(p_raw_perm = map2_dbl(vec_Sal, vec_DMSO,
                               ~ perm_vec(.x, .y, B = B, alternative = ALT)),
         p_holm = p.adjust(p_raw_perm, method = "holm"),
         p_fdr  = p.adjust(p_raw_perm, method = "BH")) %>%
  select(-vec_Sal, -vec_DMSO) %>%
  arrange(genotype, interval_num)

# Main Table S1: round-level weighted AUC sensitivity
w_star <- agg_int %>%
  group_by(genotype, interval_num) %>%
  summarise(w_star = sum(w, na.rm = TRUE), .groups = "drop")

auc_round <- agg_int %>%
  left_join(w_star, by = c("genotype", "interval_num")) %>%
  group_by(genotype, round_id, Condition) %>%
  summarise(wAUC = sum(prop * w_star, na.rm = TRUE) / sum(w_star, na.rm = TRUE),
            .groups = "drop")

sal_wauc_table <- auc_round %>%
  group_by(genotype) %>%
  summarise(stat_weighted_AUC = mean(wAUC[Condition == "Sal"]) - mean(wAUC[Condition == "DMSO"]),
            p_raw_perm = as.numeric(perm_panel_fast(wAUC, Condition, B = B, alt = ALT)["p"]),
            .groups = "drop") %>%
  mutate(p_holm = p.adjust(p_raw_perm, method = "holm"),
         p_fdr  = p.adjust(p_raw_perm, method = "BH")) %>%
  arrange(genotype)

sal_wauc_table

write_xlsx(list("Table S1 - weighted AUC" = sal_wauc_table,
                "Optional - interval-level analysis" = sal_interval_table),
           path = "output/tables/salbutamol_round_sensitivity.xlsx")