## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 150,
  fig.width = 6,
  fig.height = 4.5
  )

## ----setup, message=FALSE-----------------------------------------------------
library(smdi)
library(dplyr)
library(tibble)
library(gtsummary)
library(gt)
library(survival)
library(simsurv)
library(survminer)
library(usethis)
library(mice)
library(cardx)

# some global simulation parameters
seed_value <- 42
n <- 2500

## ----covar_generation---------------------------------------------------------
set.seed(seed_value)

# start with basic dataframe, covariates and their association with exposure
sim_covar <- tibble(
  exposure = rbinom(n = n, size = 1, prob = 0.4),
  age_num = rnorm(n, mean = 64 - 7.5*exposure, sd = 13.7),
  female_cat = rbinom(n, size = 1, prob = 0.39 - 0.05*exposure),
  ecog_cat = rbinom(n, size = 1, prob = 0.63 - 0.04*exposure),
  smoking_cat = rbinom(n, size = 1, prob = 0.45 + 0.1*exposure), 
  physical_cat = rbinom(n, size = 1, prob = 0.35 + 0.02*exposure),
  egfr_cat = rbinom(n, size = 1, prob = 0.20 + 0.07*exposure), 
  alk_cat = rbinom(n, size = 1, prob = 0.03),
  pdl1_num = rnorm(n, mean = 40 + 10*exposure, sd = 10.5),
  histology_cat = rbinom(n, size = 1, prob = 0.2),
  ses_cat = sample(x = c("1_low", "2_middle", "3_high"), size = n, replace = TRUE, prob = c(0.2 , 0.4, 0.4)),
  copd_cat = rbinom(n, size = 1, prob = 0.3 + 0.5*smoking_cat)
  ) %>%  
  # bring data in right format
  mutate(across(ends_with("num"), as.numeric)) %>% 
  mutate(across(ends_with("num"), function(x) round(x, digits = 2)))

## ----distributions_covars-----------------------------------------------------
sim_covar %>% 
  tbl_summary(by = "exposure") %>% 
  add_difference()

## ----pr(exposure)_tbl---------------------------------------------------------
exposure_form <- as.formula(paste("exposure ~ ", paste(colnames(sim_covar %>% select(-exposure)), collapse = " + ")))

exposure_fit <- glm(
  exposure_form,
  data = sim_covar,
  family = "binomial"
  )

exposure_fit %>% 
  tbl_regression(exponentiate = T)

## ----pr_treatment_assignment, fig.cap="Treatment assignment probabilities."----
# compute propensity score
exposure_plot <- sim_covar %>% 
  mutate(ps = fitted(exposure_fit))

# plot density
exposure_plot %>% 
  ggplot(aes(x = ps, fill = factor(exposure))) +
  geom_density(alpha = .5) +
  theme_bw() +
  labs(
    x = "Pr(exposure)",
    y = "Density",
    fill = "Exposed"
  )

## ----betas_outcome_generation-------------------------------------------------
betas_os <- c(
  exposure = log(1),
  age_num = log(1.05),
  female_cat = log(0.94),
  ecog_cat = log(1.25),
  smoking_cat = log(1.3),
  physical_cat = log(0.79),
  egfr_cat = log(0.5),
  alk_cat = log(0.91),
  pdl1_num = log(0.98),
  histology_cat = log(1.15)
  )

betas_os %>% 
  as.data.frame() %>% 
  transmute(logHR = round(`.`, 2)) %>% 
  rownames_to_column(var = "Covariate") %>% 
  mutate(HR = round(exp(logHR), 2)) %>% 
  gt()

## ----outcome_generation, message=FALSE----------------------------------------
set.seed(seed_value)

sim_df <- sim_covar %>% bind_cols(
  simsurv(
    dist = "exponential",
    lambdas = 0.05,
    betas = betas_os,
    x = sim_covar,
    maxt = 5 
    )
  ) %>% 
  select(-id)

## ----km_estimates-------------------------------------------------------------
km_overall <- survfit(Surv(eventtime, status) ~ 1, data = sim_df)
km_exposure <- survfit(Surv(eventtime, status) ~ exposure, data = sim_df)

tbl_survfit(
  list(km_overall, km_exposure),
  times = c(1, 5),
  label_header = "**{time} Years**"
  )

## -----------------------------------------------------------------------------
km_exposure <- survfit(Surv(eventtime, status) ~ exposure, data = sim_df)

ggsurvplot(
  km_exposure, 
  data = sim_df,
  conf.int = TRUE,
  surv.median.line = "hv",
  palette = "jco",
  xlab = "Time [Years]",
  legend.labs = c("Comparator", "Exposure of interest")
  )

## ----Cox_estimates------------------------------------------------------------
cox_lhs <- "survival::Surv(eventtime, status)"
cox_rhs <- paste(colnames(sim_covar), collapse = " + ")
cox_form = as.formula(paste(cox_lhs, "~ exposure +", cox_rhs))
  
cox_fit <- coxph(cox_form, data = sim_df)

cox_fit %>% 
  tbl_regression(exponentiate = T)

## ----export_complete_data, message=FALSE--------------------------------------
smdi_data_complete <- sim_df
use_data(smdi_data_complete, overwrite = TRUE)

## -----------------------------------------------------------------------------
# prepare a placeholder df for missing simulation
# we do not consider ses_cat
tmp <- smdi_data_complete %>% 
  select(-c(ses_cat))

# determine missingness pattern template
miss_pattern <- rep(1, ncol(tmp))

## ----mcar---------------------------------------------------------------------
# specify missingness pattern
# (0 = set to missing, 1 = remains complete)
mcar_col <- which(colnames(tmp)=="ecog_cat")
miss_pattern_mcar <- replace(miss_pattern, mcar_col, 0)

miss_prop_mcar <- .35

set.seed(42)
smdi_data_mcar <- ampute(
  data = tmp,
  prop = miss_prop_mcar,
  mech = "MCAR",
  patterns = miss_pattern_mcar,
  bycases = TRUE
  )$amp %>% 
  select(ecog_cat)

smdi_data_mcar %>% 
  select(ecog_cat) %>% 
  mutate(ecog_cat = forcats::fct_na_value_to_level(factor(ecog_cat), level = "missing")) %>% 
  tbl_summary()

## ----mar----------------------------------------------------------------------
# specify missingness pattern
# (0 = set to missing, 1 = remains complete)
mar_col <- which(colnames(tmp)=="egfr_cat")
miss_pattern_mar <- replace(miss_pattern, mar_col, 0)

# weights to compute missingness probability 
# by assigning a non-zero value
miss_weights_mar <- rep(1, ncol(tmp))
miss_weights_mar <- replace(miss_weights_mar, mar_col, 0)

miss_prop_mar <- .4

set.seed(42)
smdi_data_mar <- ampute(
  data = tmp,
  prop = miss_prop_mar,
  mech = "MAR",
  patterns = miss_pattern_mar,
  weights = miss_weights_mar,
  bycases = TRUE,
  type = "RIGHT"
  )$amp

smdi_data_mar %>% 
  select(egfr_cat) %>% 
  mutate(egfr_cat = forcats::fct_na_value_to_level(factor(egfr_cat), level = "missing")) %>%
  tbl_summary()

## ----create_mnar_v------------------------------------------------------------
# determine missingness pattern
mnar_v_col <- which(colnames(tmp)=="pdl1_num")
miss_pattern <- rep(1, ncol(tmp))
miss_pattern_mnar_v <- replace(miss_pattern, mnar_v_col, 0)

# weights to compute missingness probability 
# by assigning a non-zero value
# MNAR_v: covariate itself is only predictor
miss_weights_mnar_v <- rep(0, ncol(tmp))
miss_weights_mnar_v <- replace(miss_weights_mnar_v, mnar_v_col, 1)

miss_prop_mnar_v <- .2

set.seed(42)
smdi_data_mnar_v <- ampute(
  data = tmp,
  prop = miss_prop_mnar_v,
  mech = "MNAR",
  patterns = miss_pattern_mnar_v,
  weights = miss_weights_mnar_v,
  bycases = TRUE,
  type = "LEFT"
  )$amp

smdi_data_mnar_v %>% 
  select(pdl1_num) %>% 
  tbl_summary()

## -----------------------------------------------------------------------------
smdi_data <- smdi_data_complete %>% 
  select(-c(ecog_cat, egfr_cat, pdl1_num)) %>% 
  bind_cols(ecog_cat = smdi_data_mcar$ecog_cat, egfr_cat = smdi_data_mar$egfr_cat, pdl1_num = smdi_data_mnar_v$pdl1_num) %>% 
  mutate(across(ends_with("cat"), as.factor))

## ----distributions_covars_final-----------------------------------------------
smdi_data %>% 
  tbl_summary(by = "exposure") %>% 
  add_overall() %>% 
  add_difference()

## ----export_missing_data, message=FALSE---------------------------------------
use_data(smdi_data, overwrite = TRUE)

