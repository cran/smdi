## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  dpi = 150,
  fig.width = 6,
  fig.height = 4.5
  )

## ----setup, message=FALSE-----------------------------------------------------
# basic setup
library(smdi)
library(ggplot2)
library(survival)
library(gt)
library(broom)
library(fastDummies)
library(survival)
library(mice)
library(dplyr)
library(gridExtra)
library(tibble)
library(forcats)

## -----------------------------------------------------------------------------
# load complete dataset
smdi_data_complete <- smdi_data_complete %>% 
  dummy_columns(
    select_columns = "ses_cat",
    remove_most_frequent_dummy = TRUE,
    remove_selected_columns = TRUE
    )

# determine missingness pattern
age_col <- which(colnames(smdi_data_complete)=="age_num")
miss_pattern <- rep(1, ncol(smdi_data_complete))
miss_pattern_age <- replace(miss_pattern, age_col, 0)

# weights to compute missingness probability
# covariate itself is only predictor
miss_weights_mnar_v <- rep(0, ncol(smdi_data_complete))
miss_weights_mnar_v <- replace(miss_weights_mnar_v, age_col, 1)

miss_prop_age <- .55

set.seed(42)
smdi_data_mnar_v <- ampute(
  data = smdi_data_complete,
  prop = miss_prop_age,
  mech = "MNAR",
  patterns = miss_pattern_age,
  weights = miss_weights_mnar_v,
  bycases = TRUE,
  type = "LEFT"
  )$amp

## ----warning=FALSE------------------------------------------------------------
# plot
bind_rows(
  smdi_data_complete %>% select(age_num) %>% mutate(dataset = "Complete"),
  smdi_data_mnar_v %>% select(age_num) %>% mutate(dataset = "MNAR(value)")
  ) %>% 
  ggplot(aes(x = dataset, y = age_num)) +
  geom_violin() +
  geom_boxplot(width = 0.3) +
  stat_summary(
    fun = "mean",
    geom = "pointrange",
    color = "red"
      ) +
  labs(
    y = "Age [years]",
    x = "Cohort"
  ) +
  theme_bw()

## -----------------------------------------------------------------------------
smdi_diagnose(
  data = smdi_data_mnar_v,
  covar = "age_num",
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  ) %>% 
  smdi_style_gt()

## -----------------------------------------------------------------------------
# outcome model (see data generation script)
cox_lhs <- "survival::Surv(eventtime, status)"
covariates <- smdi_data_complete %>% 
  select(-c(exposure, eventtime, status)) %>% 
  names()

cox_rhs <- paste(covariates, collapse = " + ")
cox_form <- as.formula(paste(cox_lhs, "~ exposure +", cox_rhs))

cox_form  

## -----------------------------------------------------------------------------
# true outcome model
cox_fit_true <- coxph(cox_form, data = smdi_data_complete) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "exposure") %>% 
  select(term, estimate, conf.low, conf.high, std.error) %>% 
  mutate(analysis = "True estimate")
 
# complete case analysis
cox_fit_cc <- coxph(cox_form, data = smdi_data_mnar_v) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "exposure") %>% 
  select(term, estimate, conf.low, conf.high, std.error) %>% 
  mutate(analysis = "Complete case analysis")

# Multiple imputation (predictive mean matching)
cox_fit_imp <- mice(
  data = smdi_data_mnar_v,
  seed = 42, 
  print = FALSE
  ) %>%
  with(
    expr = coxph(formula(paste(format(cox_form), collapse = "")))
    ) %>% 
  pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE) %>% 
  filter(term == "exposure") %>% 
  select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, std.error) %>% 
  mutate(analysis = "Multiple imputation")

forest <- bind_rows(cox_fit_true, cox_fit_cc, cox_fit_imp) %>% 
  mutate(analysis = factor(analysis, levels = c("True estimate", "Complete case analysis", "Multiple imputation"))) %>% 
  ggplot(aes(y = fct_rev(analysis))) +
  geom_point(aes(x = estimate), shape = 15, size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = cox_fit_true[["estimate"]], linetype = "dashed") +
  labs(x = "Hazard ratio (95% CI)", y= "") +
  theme_bw()

table <- bind_rows(cox_fit_true, cox_fit_cc, cox_fit_imp) %>% 
  select(-term) %>% 
  relocate(analysis, .before = estimate) %>% 
  mutate(across(where(is.numeric), ~round(.x, 2)))

grid.arrange(tableGrob(table, rows = NULL), forest)

## ----echo=FALSE, eval=FALSE---------------------------------------------------
# # compute the conditional mean difference in age
# # between dataset with complete observations
# # and dataset with partially observed age covariate
# lm_form <- as.formula(paste("age_num ~ complete_dataset + exposure + ", paste(covariates[-1], collapse = "+")))
# 
# data_combined <- rbind(
#   smdi_data_complete %>% mutate(complete_dataset = 1),
#   smdi_data_mnar_v %>% mutate(complete_dataset = 0)
#   )
# 
# cond_mean_diff <- lm(lm_form, data = data_combined)$coefficients[["complete_dataset"]]
# cond_mean_diff

## ----fig.width=6--------------------------------------------------------------
# initialize method vector
method_vector <- rep("", ncol(smdi_data_mnar_v))

# columns for narfcs imputation with sensitivity parameter
mnar_imp_method <- which(colnames(smdi_data_mnar_v) == "age_num")

# update method vector
method_vector <- replace(x = method_vector, list = c(mnar_imp_method), values = c("mnar.norm"))

# modeled over a range of deltas
narfcs_modeled <- function(i){
  
  # mnar model specification ('i' is delta parameter)
  mnar_blot <- list(age_num = list(ums = paste(i)))

  narfcs_imp <- mice(
    data = smdi_data_mnar_v,
    method = method_vector,
    blots = mnar_blot,
    seed = 42, 
    print = FALSE
    ) %>% 
    with(
      expr = coxph(formula(paste(format(cox_form), collapse = "")))
      ) %>% 
    pool() %>% 
    summary(conf.int = TRUE, exponentiate = TRUE) %>% 
    filter(term == "exposure") %>% 
    select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, std.error) %>% 
    mutate(delta = i)
  
} 

narfcs_range <- lapply(
  X = seq(-25, 25, 1),
  FUN = narfcs_modeled
  )

narfcs_range_df <- do.call(rbind, narfcs_range)

reference_lines <- tibble(
  yintercept = c(cox_fit_true[[2]]),
  Reference = c("TRUE HR"),
  color = c("darkgreen")
  )

narfcs_range_df %>% 
  ggplot(aes(x = delta, y = estimate)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(x = "Sensitivity parameter δ", y = "Hazard ratio (95% CI)") +
  scale_x_continuous(breaks = seq(-25, 25, 2), limits = c(-25, 25)) +
  scale_y_continuous(breaks = seq(0.6, 1.2, 0.1), limits = c(0.6, 1.2)) +
  geom_hline(aes(yintercept = yintercept, color = Reference), reference_lines) +
  scale_colour_manual(values = reference_lines$color) +
  theme_bw() +
  theme(legend.position="top")

## -----------------------------------------------------------------------------
smdi_data %>% 
  smdi_summarize()

## -----------------------------------------------------------------------------
# we one hot encode the `ses_cat` variable again 
# in the smdi_data dataset
smdi_data <- smdi_data %>% 
  dummy_columns(
    select_columns = "ses_cat",
    remove_most_frequent_dummy = TRUE,
    remove_selected_columns = TRUE
    )

# initialize method vector
method_vector <- rep("", ncol(smdi_data))

# specify columns for narfcs and 'normal' imputation
pdl1_mnar_col <- which(colnames(smdi_data) == "pdl1_num")
ecog_mar_col <- which(colnames(smdi_data) == "ecog_cat")
egfr_mnar_col <- which(colnames(smdi_data) == "egfr_cat")

# update method vector
method_vector <- replace(
  x = method_vector, 
  list = c(pdl1_mnar_col, ecog_mar_col, egfr_mnar_col), 
  values = c("mnar.norm", "logreg", "logreg")
  )

# modeled over a range of deltas
narfcs_modeled <- function(i){
  
  # mnar model specification for pdl1_num ('i' is delta parameter)
  mnar_blot <- list(pdl1_num = list(ums = paste(i)))

  narfcs_imp <- mice(
    data = smdi_data,
    method = method_vector,
    blots = mnar_blot,
    seed = 42, 
    print = FALSE
    ) %>%
    with(
      expr = coxph(formula(paste(format(cox_form), collapse = "")))
      ) %>%
    pool() %>%
    summary(conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "exposure") %>%
    select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, std.error) %>%
    mutate(delta = i)
  
} 

narfcs_range <- lapply(
  X = seq(-25, 25, 1),
  FUN = narfcs_modeled
  )

narfcs_range_df <- do.call(rbind, narfcs_range)

reference_lines <- tibble(
  yintercept = c(cox_fit_true[[2]]),
  Reference = c("TRUE HR"),
  color = c("darkgreen")
  )

narfcs_range_df %>% 
  ggplot(aes(x = delta, y = estimate)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(x = "Sensitivity parameter δ", y = "Hazard ratio (95% CI)") +
  scale_x_continuous(breaks = seq(-25, 25, 2), limits = c(-25, 25)) +
  scale_y_continuous(breaks = seq(0.5, 1.2, 0.1), limits = c(0.5, 1.2)) +
  geom_hline(aes(yintercept = yintercept, color = Reference), reference_lines) +
  scale_colour_manual(values = reference_lines$color) +
  theme_bw() +
  theme(legend.position="top")

