#' Takes an object of class smdi and styles it to a publication-ready gt table
#'
#' `r lifecycle::badge("experimental")`
#' @description
#' This function takes either an object of class smdi or data.frame or tibble as
#' input and styles it to a publication-ready table based on the gt package.
#' The output is of class gt and can take further gt-based arguments for
#' customization.
#'
#' @param smdi_object object of class "smdi" or data.frame/tibble
#' @param include_little can be logical (TRUE/FALSE) for displaying Little's p-value that is part of an "smdi" object or a separate object of class "little"
#' @param font_size integer to determine table font size
#' @param tbl_width integer to determine table width
#'
#'
#' @return returns a formatted gt table object
#'
#' @importFrom magrittr '%>%'
#' @importFrom gt gt tab_footnote html cells_column_labels cols_label tab_options px
#' @importFrom glue glue
#'
#' @seealso
#' \code{\link[gt]{gt}}
#'
#' @export
#'
#' @examples
#'library(smdi)
#'library(dplyr)
#'
#' smdi_diagnose(
#'   data = smdi_data,
#'   covar = "egfr_cat",
#'   model = "cox",
#'   form_lhs = "Surv(eventtime, status)"
#'   ) %>%
#' smdi_style_gt()
#'

smdi_style_gt <- function(smdi_object = NULL,
                          include_little = TRUE,
                          font_size = 13,
                          tbl_width = 800
                          ){


  asmd_median_min_max <- hotteling_p <- rf_auc <- estimate_univariate <- estimate_adjusted <- NULL

  # check if smdi object or table
  if(inherits(smdi_object, "smdi")){

    smdi_table <- smdi_object$smdi_tbl

  }else if(any(class(smdi_object) %in% c("data.frame", "tibble", "tbl_df", "tbl"))){

    smdi_table <- smdi_object

  }else{

    stop("<smdi_object> is not of type smdi, data.frame or tibble")

  }

  # little checks
  if(isTRUE(include_little) & !inherits(smdi_object, "smdi")){

    warning("If include_little = TRUE, <smdi_object> needs to be of class 'smdi' (not dataframe or tibble). No p-value for Little displayed.")

  }

  # little (if specified)
  if(isTRUE(include_little) & inherits(smdi_object, "smdi")){

    # in smdi object, the p-value is already formatted
    little_foot <- glue::glue("{stringr::str_replace(smdi_object$p_little, '_', ' ')}, ")

    }else if(inherits(include_little, "little")){

    # same formatting as in smdi_diagnose():
    p_little_value <- ifelse(include_little$p.value < 0.001, '<.001', formatC(include_little$p.value, format = 'f', digits = 3))
    little_foot <- glue::glue("p little: {p_little_value}, ")

    }else{

    little_foot <- ""

    }

  # general abbrevations
  foot_abbr <- "ASMD = Median absolute standardized mean difference across all covariates, AUC = Area under the curve, beta = beta coefficient, CI = Confidence interval, max = Maximum, min = Minimum"

  smdi_gt <- smdi_table %>%

    # make into a gt table
    gt::gt() %>%

    gt::tab_footnote(
      footnote = gt::md(glue::glue("{little_foot} Abbreviations: {foot_abbr}"))
      ) %>%

    gt::cols_label(
      covariate = "Covariate",
      asmd_median_min_max= "ASMD (min/max)",
      hotteling_p = gt::md("p Hotelling"),
      rf_auc = "AUC",
      estimate_univariate = gt::md("beta univariate (95% CI)"),
      estimate_adjusted = gt::md("beta (95% CI)")
      ) %>%

    # add footnotes describing three group diagnostics
    gt::tab_footnote(
      footnote = "Group 1 diagnostic: Differences in patient characteristics between patients with and without covariate",
      locations = gt::cells_column_labels(
        columns = c(asmd_median_min_max, hotteling_p)
        )
      ) %>%

    gt::tab_footnote(
      footnote = "Group 2 diagnostic: Ability to predict missingness",
      locations = gt::cells_column_labels(
        columns = rf_auc
        )
      ) %>%

    gt::tab_footnote(
      footnote = "Group 3 diagnostic: Assessment if missingness is associated with the outcome (univariate, adjusted)",
      locations = gt::cells_column_labels(
        columns = c(estimate_univariate, estimate_adjusted)
        )
      ) %>%

    gt::tab_options(
      table.width = gt::px(tbl_width),
      data_row.padding = gt::px(3),
      table.font.size = font_size
      ) %>%

    gt::cols_align(align = "left")


  return(smdi_gt)

}

