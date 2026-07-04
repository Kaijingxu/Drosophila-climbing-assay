library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(DHARMa)
library(scales)
library(patchwork)
library(binom)
library(performance)

# Import and clean Day 21 vehicle-control data
df <- readxl::read_excel("data/day21_vehicle_only.xlsx") %>%
  tidyr::fill(`Fly No.`, Genotype, Replicates, Round, Tray, .direction = "down") %>%
  rename(fly_id   = `Fly No.`,
         genotype = Genotype,
         trial_num = Triplicate,
         replicate = Replicates,
         round_id   = Round,
         tray_local = Tray) %>%
  mutate(fly_id = factor(fly_id),
         genotype = factor(genotype,
                           levels = c("+/+", "s64w/+", "s64w/s64w")),
         trial_num = factor(trial_num, levels = c(1, 2, 3), ordered = TRUE),
         round_id = factor(round_id),
         tray_local = factor(tray_local),
         rep_chr = str_to_lower(as.character(replicate)),
         rep_chr = str_replace(rep_chr, "^rep\\s*", "r"),
         rep_num = readr::parse_number(rep_chr),
         replicate = factor(paste0("R", rep_num),
                            levels = paste0("R", sort(unique(rep_num))),
                            ordered = TRUE))

# Convert interval scores to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("Interval"),
               names_to = "interval",
               names_prefix = "Interval ",
               values_to = "success",
               names_transform = list(interval = readr::parse_number)) %>%
  mutate(interval = factor(interval, levels = 1:4, ordered = TRUE),
         interval_num = as.integer(interval),
         success = as.integer(success)) %>%
  filter(!is.na(success),
         success %in% c(0L, 1L),
         !is.na(genotype),
         !is.na(replicate))

# R1 GLMM example: Day 21 vehicle-control replicate
d_R1 <- df_long %>% filter(replicate == "R1") %>% droplevels()
m_R1 <- glmer(success ~ genotype * poly(interval_num,2) + trial_num + 
                (1 | fly_id) + 
                (1 | round_id) + 
                (1 | round_id:tray_local),
              data = d_R1, 
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5)))
summary(m_R1)

# Pearson overdispersion check
overdisp_R1 <- sum(resid(m_R1, type = "pearson")^2) / df.residual(m_R1)
cat("Pearson overdispersion ratio for R1 =", overdisp_R1, "\n")

# GLMM-predicted response probabilities by genotype and interval
emm_R1 <- emmeans(m_R1, ~ genotype | interval_num, 
               at = list(interval_num = c(1,2,3,4)),
               weights = "proportional")

plot_data_R1 <- summary(emm_R1, type = "response", infer = TRUE) %>%
  as.data.frame() %>%
  mutate(interval_num = as.integer(interval_num),
         genotype = fct_relevel(genotype, "+/+", "s64w/+", "s64w/s64w"),
         replicate = "R1")

# Contrast: homozygous mutant versus pooled wild type and heterozygote
emm_R1_resp <- regrid(emm_R1, transform = "response")

genotype_levels_R1 <- levels(emm_R1_resp)$genotype

mut_vs_pooled <- c("+/+" = -0.5,
                   "s64w/+" = -0.5,
                   "s64w/s64w" = 1)

mut_vs_pooled_R1 <- setNames(mut_vs_pooled[genotype_levels_R1],
                             genotype_levels_R1)

cmp_R1 <- contrast(emm_R1_resp,
                   method = list("s64w/s64w - pooled(+/+ and s64w/+)" = mut_vs_pooled_R1),
                   by = "interval_num",
                   adjust = "none")

out_R1 <- summary(cmp_R1, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(interval_num = as.integer(interval_num),
         RD_pct = estimate * 100,
         LCL_pct = asymp.LCL * 100,
         UCL_pct = asymp.UCL * 100,
         p_holm = p.adjust(p.value, method = "holm"),
         label = if_else(p_holm < 1e-3,
                         sprintf("p = %.1e", p_holm),
                         paste0("p = ", formatC(p_holm, digits = 3, format = "f")))) %>%
  select(interval_num, contrast,
         estimate, SE, asymp.LCL, asymp.UCL,
         RD_pct, LCL_pct, UCL_pct,
         z.ratio, p.value, p_holm, label)

# fly number by genotype 
fly_n_R1 <- d_R1 %>%
  distinct(genotype, fly_id) %>%
  count(genotype, name = "n_fly") %>%
  pivot_wider(names_from = genotype,
              values_from = n_fly,
              names_prefix = "n_")

# Marginal predictions for calibration
pred_R1 <- d_R1 %>%
  mutate(pred = predict(m_R1, newdata = d_R1, type = "response", re.form = NA))

# Binned calibration with Wilson confidence intervals
nbins <- 8
bin_breaks_R1 <- unique(quantile(pred_R1$pred,
                                 probs = seq(0, 1, length.out = nbins + 1),
                                 type = 8,
                                 na.rm = TRUE))

if (length(bin_breaks_R1) < 3) {
  pred_range_R1 <- range(pred_R1$pred, na.rm = TRUE)
  bin_breaks_R1 <- seq(pred_range_R1[1],
                       pred_range_R1[2],
                       length.out = min(nbins + 1, 3))}

pred_R1 <- pred_R1 %>%
  mutate(pred_bin = cut(pred,
                        breaks = bin_breaks_R1,
                        include.lowest = TRUE,
                        right = FALSE))

calib_R1 <- pred_R1 %>%
  group_by(pred_bin) %>%
  summarise(n = n(),
            pred_mean = mean(pred, na.rm = TRUE),
            obs_rate = mean(success, na.rm = TRUE),
            successes = sum(success, na.rm = TRUE),
            .groups = "drop")

calib_ci_R1 <- binom::binom.confint(x = calib_R1$successes,
                                    n = calib_R1$n,
                                    methods = "wilson")

calib_R1 <- bind_cols(calib_R1,
                      calib_ci_R1[, c("lower", "upper")])

# Calibration slope, intercept and Brier score
eps <- 1e-6
pred_R1 <- pred_R1 %>%
  mutate(pred_clip = pmin(pmax(pred, eps), 1 - eps),
         logit_pred = qlogis(pred_clip))

fit_cal_R1 <- glm(success ~ logit_pred,
                  data = pred_R1,
                  family = binomial(link = "logit"))

cal_summary_R1 <- tibble(replicate = "R1",
                         calibration_slope = coef(fit_cal_R1)[["logit_pred"]],
                         calibration_intercept = coef(fit_cal_R1)[["(Intercept)"]],
                         brier = mean((pred_R1$success - pred_R1$pred)^2, na.rm = TRUE),
                         overdispersion_ratio = overdisp_R1)

cal_label_R1 <- sprintf("slope = %.2f, intercept = %.2f, Brier = %.3f",
                        cal_summary_R1$calibration_slope,
                        cal_summary_R1$calibration_intercept,
                        cal_summary_R1$brier)

# R1-R5 running
run_vehicle_replicate <- function(rep_id, data = df_long) {
  d_rep <- data %>%
    filter(replicate == rep_id) %>%
    droplevels()
  
  m_rep <- glmer(
    success ~ genotype * poly(interval_num, 2) + trial_num +
      (1 | fly_id) +
      (1 | round_id) +
      (1 | round_id:tray_local),
    data = d_rep,
    family = binomial(link = "logit"),
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  overdisp <- sum(resid(m_rep, type = "pearson")^2) / df.residual(m_rep)
  
  emm_rep <- emmeans(
    m_rep,
    ~ genotype | interval_num,
    at = list(interval_num = 1:4),
    weights = "proportional"
  )
  
  plot_data <- summary(emm_rep, type = "response", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(
      interval_num = as.integer(interval_num),
      genotype = fct_relevel(genotype, "+/+", "s64w/+", "s64w/s64w"),
      replicate = rep_id
    )
  
  emm_rep_resp <- regrid(emm_rep, transform = "response")
  
  genotype_levels <- levels(emm_rep_resp)$genotype
  
  mut_vs_pooled <- c(
    "+/+" = -0.5,
    "s64w/+" = -0.5,
    "s64w/s64w" = 1
  )
  
  mut_vs_pooled <- setNames(mut_vs_pooled[genotype_levels], genotype_levels)
  
  cmp <- contrast(
    emm_rep_resp,
    method = list("s64w/s64w - pooled(+/+ and s64w/+)" = mut_vs_pooled),
    by = "interval_num",
    adjust = "none"
  )
  
  contrast_data <- summary(cmp, infer = TRUE) %>%
    as.data.frame() %>%
    mutate(
      interval_num = as.integer(interval_num),
      replicate = rep_id,
      RD_pct = estimate * 100,
      LCL_pct = asymp.LCL * 100,
      UCL_pct = asymp.UCL * 100,
      p_holm = p.adjust(p.value, method = "holm"),
      label = if_else(
        p_holm < 1e-3,
        sprintf("p = %.1e", p_holm),
        paste0("p = ", formatC(p_holm, digits = 3, format = "f"))
      )
    )
  
  pred_data <- d_rep %>%
    mutate(
      pred = predict(
        m_rep,
        newdata = d_rep,
        type = "response",
        re.form = NA
      )
    )
  
  eps <- 1e-6
  
  pred_data <- pred_data %>%
    mutate(
      pred_clip = pmin(pmax(pred, eps), 1 - eps),
      logit_pred = qlogis(pred_clip)
    )
  
  fit_cal <- glm(
    success ~ logit_pred,
    data = pred_data,
    family = binomial(link = "logit")
  )
  
  cal_summary <- tibble(
    replicate = rep_id,
    overdispersion_ratio = overdisp,
    calibration_slope = coef(fit_cal)[["logit_pred"]],
    calibration_intercept = coef(fit_cal)[["(Intercept)"]],
    brier = mean((pred_data$success - pred_data$pred)^2, na.rm = TRUE)
  )
  
  list(
    model = m_rep,
    emmeans = emm_rep,
    plot_data = plot_data,
    contrast_data = contrast_data,
    prediction_data = pred_data,
    calibration_data = calib_R1,
    calibration_summary = cal_summary
  )
}

vehicle_results <- levels(df_long$replicate) %>%
  set_names() %>%
  map(run_vehicle_replicate, data = df_long)

vehicle_plot_data <- map_dfr(vehicle_results, "plot_data")
vehicle_contrasts <- map_dfr(vehicle_results, "contrast_data")
vehicle_calibration <- map_dfr(vehicle_results, "calibration_summary")
vehicle_models <- map(vehicle_results, "model")

# Figure 3: per-replicate vehicle-control GLMM predictions
replicate_levels <- paste0("R", 1:5)

plot_data_all <- vehicle_plot_data %>%
  mutate(genotype = fct_relevel(genotype, "+/+", "s64w/+", "s64w/s64w"),
         interval_num = as.integer(interval_num),
         replicate = factor(replicate, levels = replicate_levels, ordered = TRUE))

contrast_all <- vehicle_contrasts %>%
  mutate(interval_num = as.integer(interval_num),
         replicate = factor(replicate, levels = replicate_levels, ordered = TRUE))

vehicle_map <- c(R1 = "Water",
                 R2 = "Water",
                 R3 = "0.5% DMSO",
                 R4 = "0.016% DMSO",
                 R5 = "0.016% DMSO")

replicate_labels <- df_long %>%
  mutate(vehicle = vehicle_map[as.character(replicate)]) %>%
  distinct(replicate, fly_id, vehicle) %>%
  count(replicate, vehicle, name = "n") %>%
  mutate(replicate = factor(replicate, levels = replicate_levels, ordered = TRUE),
         label = paste0(as.character(replicate), " (", vehicle, "; n = ", n, ")"))

replicate_label_vector <- setNames(replicate_labels$label,
                                   as.character(replicate_labels$replicate))

label_positions <- plot_data_all %>%
  left_join(contrast_all %>%
              transmute(replicate,
                        interval_num,
                        label),
            by = c("replicate", "interval_num")) %>%
  mutate(headroom = 1 - asymp.UCL,
         y = case_when(headroom >= 0.07 ~ pmin(0.98, asymp.UCL + 0.05),
                       prob >= 0.07     ~ pmax(0.02, prob - 0.05),
                       TRUE             ~ pmin(0.98, prob + 0.05)),
         label_show = if_else(genotype == "s64w/s64w", label, ""))

pd <- position_dodge(width = 0.25)

p_Figure_three <- ggplot(plot_data_all,
                   aes(x = interval_num,
                       y = prob,
                       color = genotype,
                       shape = genotype,
                       group = genotype)) +
  geom_line(linewidth = 0.9, position = pd) +
  geom_point(size = 2.7, position = pd) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.10,
                position = pd) +
  facet_wrap(~ replicate,
             nrow = 1,
             labeller = labeller(replicate = replicate_label_vector)) +
  scale_x_continuous(breaks = 1:4,
                     labels = paste("Interval", 1:4)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0, 1.05), clip = "off") +
  scale_color_manual(values = c("+/+" = "blue",
                                "s64w/+" = "green",
                                "s64w/s64w" = "red")) +
  scale_shape_manual(values = c("+/+" = 16,
                                "s64w/+" = 17,
                                "s64w/s64w" = 15)) +
  labs(title = "Per-replicate GLMM-predicted success in Day 21 vehicle controls",
       x = "Interval",
       y = "Predicted probability of success",
       color = "Genotype",
       shape = "Genotype") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold", size = 16),
        legend.title = element_text(face = "bold")) +
  geom_text(data = label_positions,
            aes(x = interval_num, y = y, label = label_show, group = genotype),
            position = pd,
            inherit.aes = FALSE,
            size = 3,
            vjust = 0)

p_Figure_three


# Supplementary figure 2: apparent calibration of Day 21 vehicle-control GLMMs
calib_long <- map_dfr(vehicle_results,
                      ~ .x$calibration_data,
                      .id = "replicate") %>%
  mutate(replicate = factor(replicate, levels = paste0("R", 1:5)))

calib_labels <- vehicle_calibration %>%
  mutate(replicate = factor(replicate, levels = paste0("R", 1:5)),
         x = 0.02,
         y = 0.98,
         label = sprintf("slope = %.2f, int = %.2f, Brier = %.3f",
                         calibration_slope,
                         calibration_intercept,
                         brier))

# plot
p_cal_facet <- ggplot(calib_long, aes(x = pred_mean, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  facet_wrap(~ replicate, nrow = 2) +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = 0.02)) +
  coord_equal() +
  labs(title = "Apparent calibration of Day 21 vehicle-control GLMMs",
       x = "Predicted probability (binned mean)",
       y = "Observed success rate (Wilson CI)") +
  geom_text(data = calib_labels,
            aes(x = x, y = y, label = label),
            hjust = 0,
            vjust = 1,
            size = 3.2,
            inherit.aes = FALSE) +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(14, "pt"),
        panel.spacing.y = unit(10, "pt"),
        plot.margin = margin(5.5, 12, 5.5, 12))

p_cal_facet


# Supplementary Figure 3A-C: diagnostics for Day 21 vehicle-control GLMMs
replicate_levels <- paste0("R", 1:5)
vehicle_map <- c(R1 = "Water",
                 R2 = "Water",
                 R3 = "0.5% DMSO",
                 R4 = "0.016% DMSO",
                 R5 = "0.016% DMSO")

vehicle_models <- map(vehicle_results, "model")

# Helper: safe Nakagawa R2 extraction
r2_safe <- function(model) {
  r2_out <- try(performance::r2_nakagawa(model), silent = TRUE)
  
  if (!inherits(r2_out, "try-error")) {
    r2_marginal <- suppressWarnings(as.numeric(r2_out$R2_marginal))
    r2_conditional <- suppressWarnings(as.numeric(r2_out$R2_conditional))
    
    if (!is.na(r2_marginal)) {
      if (is.na(r2_conditional)) r2_conditional <- r2_marginal
      return(list(R2_marginal = r2_marginal, R2_conditional = r2_conditional))
    }
  }
  
  eta_fixed <- try(predict(model, type = "link", re.form = NA), silent = TRUE)
  
  if (inherits(eta_fixed, "try-error")) {
    X <- try(lme4::getME(model, "X"), silent = TRUE)
    beta <- try(lme4::fixef(model), silent = TRUE)
    
    eta_fixed <- if (!inherits(X, "try-error") && !inherits(beta, "try-error")) {
      as.numeric(X %*% beta)
    } else {
      NA_real_
    }
  }
  
  var_fixed <- stats::var(eta_fixed, na.rm = TRUE)
  
  vc <- as.data.frame(lme4::VarCorr(model))
  sd_intercepts <- vc$sdcor[vc$var1 == "(Intercept)"]
  var_random <- sum(sd_intercepts^2, na.rm = TRUE)
  var_residual <- (pi^2) / 3
  
  denominator <- var_fixed + var_random + var_residual
  
  if (!is.finite(denominator) || denominator <= 0) {
    return(list(R2_marginal = NA_real_, R2_conditional = NA_real_))
  }
  
  list(
    R2_marginal = var_fixed / denominator,
    R2_conditional = (var_fixed + var_random) / denominator
  )
}

get_model_metrics <- function(model, replicate_id) {
  overdispersion <- sum(resid(model, type = "pearson")^2) / df.residual(model)
  
  r2 <- r2_safe(model)
  
  vc <- as.data.frame(lme4::VarCorr(model))
  
  get_sd <- function(group_name) {
    x <- vc %>%
      filter(grp == group_name, var1 == "(Intercept)") %>%
      pull(sdcor)
    
    if (length(x) == 0) 0 else x
  }
  
  sd_fly <- get_sd("fly_id")
  sd_round <- get_sd("round_id")
  sd_tray <- get_sd("round_id:tray_local")
  
  var_fly <- sd_fly^2
  var_round <- sd_round^2
  var_tray <- sd_tray^2
  var_residual <- (pi^2) / 3
  denominator <- var_fly + var_round + var_tray + var_residual
  
  set.seed(123)
  sim_res <- suppressWarnings(DHARMa::simulateResiduals(model, n = 1000))
  
  tibble(
    replicate = factor(replicate_id, levels = replicate_levels, ordered = TRUE),
    vehicle = vehicle_map[[replicate_id]],
    overdispersion = overdispersion,
    R2_marginal = r2$R2_marginal,
    R2_conditional = r2$R2_conditional,
    SD_fly = sd_fly,
    SD_round = sd_round,
    SD_tray = sd_tray,
    ICC_fly = var_fly / denominator,
    ICC_round = var_round / denominator,
    ICC_tray = var_tray / denominator,
    p_uniform = suppressWarnings(DHARMa::testUniformity(sim_res)$p.value),
    p_dispersion = suppressWarnings(DHARMa::testDispersion(sim_res)$p.value),
    p_outliers = suppressWarnings(DHARMa::testOutliers(sim_res)$p.value),
    singular = lme4::isSingular(model, tol = 1e-5)
  )
}

vehicle_metrics <- imap_dfr(vehicle_models, get_model_metrics) %>%
  mutate(R2_marginal = pmax(0, pmin(1, R2_marginal)),
         R2_conditional = if_else(is.na(R2_conditional),
                                  R2_marginal,
                                  R2_conditional),
         R2_conditional = pmax(0, pmin(1, R2_conditional)))

### Supplementary Figure 3A: DHARMa tests
dharma_long <- vehicle_metrics %>%
  select(replicate, p_uniform, p_outliers) %>%
  pivot_longer(cols = -replicate,
               names_to = "test",
               values_to = "p_value") %>%
  mutate(test = recode(test,
                       p_uniform = "Uniformity",
                       p_outliers = "Outliers"),
         test = factor(test, levels = c("Uniformity", "Outliers")),
         significant = p_value < 0.05)

p_dharma <- ggplot(dharma_long, aes(x = replicate, y = test)) +
  geom_point(aes(shape = significant), size = 3, stroke = 1.1) +
  geom_text(aes(label = sprintf("%.3f", p_value)), nudge_y = 0.20, size = 4) +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "p ≥ 0.05"),
                     guide = guide_legend(title = NULL)) +
  labs(title = "A: DHARMa tests (simulation p-values)",
       x = NULL,
       y = NULL) +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"))

p_dharma

### Supplementary Figure 3B: Pearson dispersion ratios
vehicle_metrics <- vehicle_metrics %>%
  mutate(overdisp_display = round(overdispersion, 2),
         dispersion_class = case_when(overdisp_display < 0.50 ~ "<0.5 (under)",
                                      overdisp_display > 2.00 ~ ">2 (over)",
                                      TRUE                    ~ "0.5–2 (OK)"))

p_overdisp <- ggplot(vehicle_metrics, aes(x = replicate, y = overdispersion, shape = dispersion_class)) +
  annotate("rect",
           xmin = 0.5,
           xmax = length(replicate_levels) + 0.5,
           ymin = 0.5,
           ymax = 2,
           fill = "grey95") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", overdispersion)), vjust = -0.9, size = 4) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25", "0.5", "1", "2", "4"),
                limits = c(0.25, 4)) +
  scale_shape_manual(values = c("0.5–2 (OK)" = 16, "<0.5 (under)" = 4, ">2 (over)" = 4)) +
  labs(title = "B: Pearson dispersion ratios",
       x = NULL,
       y = "Ratio (ideal = 1; acceptable 0.5-2)",
       shape = NULL) +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"))

p_overdisp

### Supplementary Figure 3C: random-intercept BLUP distributions
extract_random_effects <- function(model, group_name, replicate_id) {
  random_effects <- ranef(model, condVar = TRUE)
  
  if (!(group_name %in% names(random_effects))) {
    return(tibble())
  }
  
  group_df <- random_effects[[group_name]]
  effects <- group_df[["(Intercept)"]]
  post_var <- attr(group_df, "postVar")
  
  se <- sqrt(sapply(seq_along(effects), function(i) post_var[1, 1, i]))
  
  sd_group <- as.data.frame(lme4::VarCorr(model)) %>%
    filter(grp == group_name, var1 == "(Intercept)") %>%
    pull(sdcor)
  
  tibble(replicate = replicate_id,
         group = if_else(group_name == "round_id", "Round", "Tray within round"),
         level_id = rownames(group_df),
         effect = as.numeric(effects),
         SE = se,
         lower = effect - 1.96 * SE,
         upper = effect + 1.96 * SE,
         sd_model = sd_group)}

random_effects_all <- imap_dfr(vehicle_models, function(model, replicate_id) {
  bind_rows(extract_random_effects(model, "round_id", replicate_id),
            extract_random_effects(model, "round_id:tray_local", replicate_id))}) %>%
  mutate(replicate = factor(replicate, levels = replicate_levels, ordered = TRUE),
         group = factor(group, levels = c("Round", "Tray within round")))

random_effect_labels <- random_effects_all %>%
  group_by(replicate, group) %>%
  summarise(SD = first(sd_model),
            n_units = n(),
            .groups = "drop") %>%
  mutate(label = if_else(SD < 1e-3,
                         paste0("SD ≈ 0 (boundary); n = ", n_units),
                         paste0("SD = ", sprintf("%.2f", SD), "; n = ", n_units)))

y_min <- min(random_effects_all$lower, na.rm = TRUE)
y_max <- max(random_effects_all$upper, na.rm = TRUE)
y_pad <- 0.1 * (y_max - y_min + 1e-8)
y_limits <- c(y_min - y_pad, y_max + y_pad)

label_y_positions <- random_effects_all %>%
  group_by(replicate, group) %>%
  summarise(y_position = max(upper, na.rm = TRUE) + 0.08 * diff(y_limits),
            .groups = "drop")

p_random_effects <- ggplot(random_effects_all, aes(x = replicate, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "grey95", color = "grey70") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, alpha = 0.35) +
  geom_point(position = position_jitter(width = 0.08, height = 0),
             alpha = 0.7,
             size = 1.6) +
  facet_grid(rows = vars(group), scales = "fixed") +
  coord_cartesian(ylim = y_limits) +
  labs(title = "C: Random-intercept BLUP distributions by replicate",
       x = "Replicate (vehicle-control, Day21)",
       y = "Random intercept (log-odds, BLUP)") +
  geom_text(data = random_effect_labels %>%
              left_join(label_y_positions, by = c("replicate", "group")),
            aes(x = replicate, y = y_position, label = label),
            inherit.aes = FALSE,
            size = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"))

p_random_effects

# Supplementary Figure 4A-B: GLMM vs ablation GLM in R1
log_loss <- function(y, p) {eps <- 1e-12
-mean(y * log(p + eps) + (1 - y) * log(1 - p + eps))}

brier_score <- function(y, p) {mean((p - y)^2)}

bin_calibration <- function(data, pred_col = "pred", outcome_col = "success", nbins = 8) {
  pred <- data[[pred_col]]
  outcome <- data[[outcome_col]]
  
  breaks <- unique(quantile(pred, probs = seq(0, 1, length.out = nbins + 1), type = 8, na.rm = TRUE))
  
  if (length(breaks) < 3) {pred_range <- range(pred, na.rm = TRUE)
  breaks <- seq(pred_range[1], pred_range[2], length.out = min(nbins + 1, 3))}
  
  calib <- data %>%
    mutate(pred_bin = cut(.data[[pred_col]], breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
    group_by(pred_bin) %>%
    summarise(n = n(),
              pred_mean = mean(.data[[pred_col]], na.rm = TRUE),
              obs_rate = mean(.data[[outcome_col]], na.rm = TRUE),
              successes = sum(.data[[outcome_col]], na.rm = TRUE),
              .groups = "drop")
  
  ci <- binom::binom.confint(x = calib$successes,
                             n = calib$n,
                             methods = "wilson")
  
  bind_cols(calib, ci[, c("lower", "upper")])
}

calibration_stats <- function(y, p) {eps <- 1e-6
p_clip <- pmin(pmax(p, eps), 1 - eps)
fit <- glm(y ~ qlogis(p_clip), family = binomial(link = "logit"))

  tibble(log_loss = log_loss(y, p),
         brier = brier_score(y, p),
         cal_intercept = unname(coef(fit)[1]),
         cal_slope = unname(coef(fit)[2]))}

# R1 data and models
d_R1 <- df_long %>%
  filter(replicate == "R1") %>%
  droplevels()

m_R1_glmm <- vehicle_results$R1$model

m_R1_glm <- glm(success ~ genotype * poly(interval_num, 2) + trial_num,
                data = d_R1,
                family = binomial(link = "logit"))

pred_R1_glmm <- vehicle_results$R1$prediction_data %>%
  select(success, pred) %>%
  mutate(model = "GLMM")

pred_R1_glm <- d_R1 %>%
  mutate(pred = predict(m_R1_glm, newdata = d_R1, type = "response"),
         model = "Ablation GLM") %>%
  select(success, pred, model)

calib_R1_glmm <- vehicle_results$R1$calibration_data %>%
  mutate(model = "GLMM")

calib_R1_glm <- bin_calibration(pred_R1_glm) %>%
  mutate(model = "Ablation GLM")

stats_R1_app <- bind_rows(calibration_stats(pred_R1_glmm$success, pred_R1_glmm$pred) %>% mutate(model = "GLMM"),
                          calibration_stats(pred_R1_glm$success, pred_R1_glm$pred) %>% mutate(model = "Ablation GLM")) %>%
  mutate(label = sprintf("slope = %.2f, int = %.2f, Brier = %.3f",
                         cal_slope, cal_intercept, brier))

# Supplementary Figure 4A: apparent calibration
calib_R1_app <- bind_rows(calib_R1_glmm, calib_R1_glm) %>%
  mutate(model = factor(model, levels = c("GLMM", "Ablation GLM")))

calib_R1_app <- calib_R1_app %>%
  mutate(model_label = recode(model,
                              "GLMM" = "GLMM (fixed-effects predictions; re.form = NA)",
                              "Ablation GLM" = "Ablation GLM (no random effects)"),
         model_label = factor(model_label,
                              levels = c("GLMM (fixed-effects predictions; re.form = NA)",
                                         "Ablation GLM (no random effects)")))

stats_R1_app <- stats_R1_app %>%
  mutate(model_label = recode(model,
                              "GLMM" = "GLMM (fixed-effects predictions; re.form = NA)",
                              "Ablation GLM" = "Ablation GLM (no random effects)"),
         model_label = factor(model_label,
                              levels = levels(calib_R1_app$model_label)))

stats_R1_app <- stats_R1_app %>%
  mutate(model = factor(model, levels = c("GLMM", "Ablation GLM")),
         x = 0.02,
         y = 0.98)

p_supp4A <- ggplot(calib_R1_app, aes(x = pred_mean, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  facet_wrap(~ model_label, nrow = 1) +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = 0.02)) +
  coord_equal() +
  geom_text(data = stats_R1_app,
            aes(x = x, y = y, label = label),
            hjust = 0,
            vjust = 1,
            size = 3.2,
            inherit.aes = FALSE) +
  labs(title = "A: R1 apparent calibration: GLMM vs ablation GLM",
       x = "Predicted probability (binned mean)",
       y = "Observed success rate (Wilson CI)") +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"),
        panel.grid.minor = element_blank())

p_supp4A

# Round-wise cross-validation
set.seed(123)

rounds <- sample(unique(d_R1$round_id))
K <- min(5L, length(rounds))
folds <- split(rounds, cut(seq_along(rounds), breaks = K, labels = FALSE))

fit_cv_model <- function(test_rounds, model_type) {
  train_data <- d_R1 %>% filter(!round_id %in% test_rounds)
  test_data  <- d_R1 %>% filter(round_id %in% test_rounds)
  
  if (model_type == "GLMM") {
    fit <- glmer(success ~ genotype * poly(interval_num, 2) + trial_num +
                   (1 | fly_id) +
                   (1 | round_id) +
                   (1 | round_id:tray_local),
                 data = train_data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5)))
    
    pred <- predict(fit, newdata = test_data, type = "response", re.form = NA)
  } else {
    fit <- glm(success ~ genotype * poly(interval_num, 2) + trial_num,
               data = train_data,
               family = binomial(link = "logit"))
    
    pred <- predict(fit, newdata = test_data, type = "response")}
  
  tibble(
    y = test_data$success,
    p = pred,
    model = model_type
  )
}

cv_predictions <- imap_dfr(folds, function(test_rounds, fold_id) {
  bind_rows(fit_cv_model(test_rounds, "GLMM"),
            fit_cv_model(test_rounds, "Ablation GLM")) %>%
    mutate(fold = as.integer(fold_id))
})

ece <- function(y, p, nbins = 10) {
  tibble(y = y, p = p) %>%
    mutate(bin = ntile(p, nbins)) %>%
    group_by(bin) %>%
    summarise(n = n(),
              pred_mean = mean(p),
              obs_rate = mean(y),
              .groups = "drop") %>%
    summarise(ece = sum(n * abs(obs_rate - pred_mean)) / sum(n)) %>%
    pull(ece)}

cv_metrics <- cv_predictions %>%
  group_by(model) %>%
  summarise(log_loss = log_loss(y, p),
            brier = brier_score(y, p),
            ece = ece(y, p),
            .groups = "drop") %>%
  mutate(model = factor(model, levels = c("Ablation GLM", "GLMM")))

# Supplementary Figure 4B: cross-validated metrics
p_supp4B <- cv_metrics %>%
  pivot_longer(cols = c(log_loss, brier, ece),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("log_loss", "brier", "ece"),
                         labels = c("Log loss", "Brier score", "ECE"))) %>%
  ggplot(aes(x = metric, y = value, fill = model)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.55) +
  labs(title = "B: R1 round-wise cross-validated metrics (CV by round)",
       x = NULL,
       y = NULL,
       fill = NULL) +
  theme_bw() +
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"),
        panel.grid.minor = element_blank())

p_supp4B

# Optional model summaries: Nakagawa R2 and ICC components
r2_icc_summary <- vehicle_metrics %>%
  select(replicate,
         vehicle,
         R2_marginal,
         R2_conditional,
         ICC_fly,
         ICC_round,
         ICC_tray) %>%
  mutate(across(c(R2_marginal, R2_conditional, ICC_fly, ICC_round, ICC_tray),
                ~ round(.x, 3)))
r2_icc_summary