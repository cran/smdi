## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 150,
  fig.width = 6,
  fig.height = 4.5
  )

## ----setup--------------------------------------------------------------------
library(smdi)
library(gt)
suppressPackageStartupMessages(library(dplyr))

## ---- fig.cap="Illustrating missing indicator variable generation within `smdi` functions"----
smdi_data %>% 
  smdi_na_indicator(
    drop_NA_col = FALSE # usually TRUE, but for demonstration purposes set to FALSE
    ) %>% 
  select(
    ecog_cat, ecog_cat_NA, 
    egfr_cat, egfr_cat_NA, 
    pdl1_num, pdl1_num_NA
    ) %>% 
  head() %>% 
  gt()

## -----------------------------------------------------------------------------
# we simulatea monotone missingness pattern
# following an MCAR mechanism
data_monotone <- smdi_data_complete %>% 
  mutate(
    lab1 = rnorm(nrow(smdi_data_complete), mean = 5, sd = 0.5),
    lab2 = rnorm(nrow(smdi_data_complete), mean = 10, sd = 2.25)
    )

data_monotone[3:503, "lab1"] <- NA
data_monotone[1:500, "lab2"] <- NA

## -----------------------------------------------------------------------------
smdi::gg_miss_upset(data = data_monotone)

## -----------------------------------------------------------------------------
smdi::md.pattern(data_monotone[, c("lab1", "lab2")], plot = FALSE)

## -----------------------------------------------------------------------------
diagnostics_jointly <- smdi_diagnose(
  data = data_monotone,
  covar = NULL, # NULL includes all covariates with at least one NA
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  )

## ---- fig.cap="Diagnostics of lab 1 if analyzed separately."------------------
diagnostics_jointly %>% 
  smdi_style_gt()

## ---- fig.cap="Diagnostics of lab 1 if analyzed separately."------------------
# lab 1
lab1_diagnostics <- smdi_diagnose(
  data = data_monotone %>% select(-lab2),
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  )

lab1_diagnostics %>% 
  smdi_style_gt()

## ---- fig.cap="Diagnostics of lab 2 if analyzed separately."------------------
# lab 2
lab2_diagnostics <- smdi_diagnose(
  data = data_monotone %>% select(-lab1),
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  )

lab2_diagnostics %>% 
  smdi_style_gt()

## -----------------------------------------------------------------------------
# computing a gloabl p-value for Little's test including both lab1 and lab2
little_global <- smdi_little(data = data_monotone)

# combining two individual lab smdi tables and global Little's test
smdi_style_gt(
  smdi_object = rbind(lab1_diagnostics$smdi_tbl, lab2_diagnostics$smdi_tbl), 
  include_little = little_global
  )

