## ----include = FALSE----------------------------------------------------------
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

## -----------------------------------------------------------------------------
smdi_diagnose(
  data = smdi_data,
  covar = NULL, # NULL includes all covariates with at least one NA
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  ) %>% 
  smdi_style_gt()

