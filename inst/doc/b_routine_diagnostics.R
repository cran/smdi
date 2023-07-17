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

## ---- fig-group_diagnostics, fig.cap = "Overview three group diagnostics", echo = FALSE----
knitr::include_graphics(here::here("vignettes", "smdi_diagnose_table.png"))

## ---- fig-guidance, fig.cap = "Example of how `smdi` diagnostics can be applied to give insights into the likelihood of underlying missingness structures in a real-world database study.", echo = FALSE----
knitr::include_graphics(here::here("vignettes", "smdi_examples.png"))

## -----------------------------------------------------------------------------
smdi_data %>% 
  dplyr::glimpse()

## ---- eval=FALSE--------------------------------------------------------------
#  # dataset with simulated missingness
#  ?smdi::smdi_data()
#  
#  # complete dataset
#  ?smdi::smdi_data_complete()

## -----------------------------------------------------------------------------
smdi_data %>% 
  smdi_summarize()

## -----------------------------------------------------------------------------
covars_missing <- smdi_summarize(data = smdi_data) %>% 
  pull(covariate)

smdi_data %>% 
  smdi_vis(covar = covars_missing)

## -----------------------------------------------------------------------------
smdi_data %>% 
  smdi_vis(covar = covars_missing, strata = "exposure")

## -----------------------------------------------------------------------------
smdi::gg_miss_upset(data = smdi_data)

## ---- fig.width = 8, fig.height=6, fig.width = 10-----------------------------
smdi::md.pattern(smdi_data[, c(covars_missing)], plot = FALSE)

## -----------------------------------------------------------------------------
asmd <- smdi_asmd(data = smdi_data)

## -----------------------------------------------------------------------------
asmd$egfr_cat$asmd_table1

## -----------------------------------------------------------------------------
asmd$egfr_cat$asmd_plot

## -----------------------------------------------------------------------------
asmd$egfr_cat$asmd_aggregate

## -----------------------------------------------------------------------------
summary(asmd)

## -----------------------------------------------------------------------------
h0 <- smdi_hotelling(data = smdi_data)
h0

## -----------------------------------------------------------------------------
h0$ecog_cat

## -----------------------------------------------------------------------------
h0_global <- smdi_little(data = smdi_data)
h0_global

## -----------------------------------------------------------------------------
auc <- smdi_rf(data = smdi_data)
auc$ecog_cat$rf_table
auc$ecog_cat$rf_plot

## ---- eval=FALSE--------------------------------------------------------------
#  ?smdi::smdi_outcome()

## -----------------------------------------------------------------------------
smdi_outcome(
  data = smdi_data, 
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  )

## -----------------------------------------------------------------------------
diagnostics <- smdi_diagnose(
  data = smdi_data,
  covar = NULL, # NULL includes all covariates with at least one NA
  model = "cox",
  form_lhs = "Surv(eventtime, status)"
  )

## -----------------------------------------------------------------------------
diagnostics$smdi_tbl

## -----------------------------------------------------------------------------
diagnostics$p_little

## -----------------------------------------------------------------------------
library(gt)

smdi_style_gt(diagnostics)

## ---- eval = FALSE------------------------------------------------------------
#  gtsave(
#    data = smdi_style_gt(diagnostics),
#    filename = "smdi_table.docx", # name of the final .docx file
#    path = "." # path where the file should be stored
#    )

