# Unit tests for the smdi_outcome function
test_that("smdi_outcome returns the expected output", {
  # Test case 1: Check if function throws an error when no dataframe is provided
  expect_error(smdi_outcome())
})

# n_cores on windows
test_that("Test n_cores OS dependency", {
  set.seed(42)
  data <- data.frame(x = c(NA, NA, rbinom(97, 1, 0.5), NA), y = c(rbinom(97, 1, 0.5), NA, NA, NA), z = rnorm(100))
  if(isTRUE(Sys.info()[["sysname"]]=="Windows")){
    expect_warning(smdi_outcome(data = data, model = "linear", form_lhs = "z", n_cores = 2))
  }else{
    expect_no_warning(smdi_outcome(data = data, model = "linear", form_lhs = "z", n_cores = 2))
  }
})

test_that("No LHS form provided", {
  # Test case 2: Check if function throws an error when no form_lhs is provided
  expect_error()
})

test_that("Invalid LHS provided", {
  # Test case 3: Check if function throws an error when an invalid model type is specified
  expect_error(smdi_outcome(data = iris, form_lhs = "Sepal.Length", model = "invalid_model"), "<model> either not specified or not of type glm, linear or cox")
})


# Test case for glm model with binary outcome
test_that("glm model with binary outcome works correctly", {
  set.seed(42)
  data <- data.frame(outcome = rbinom(100, 1, 0.5),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  result <- smdi_outcome(data = data, model = "glm", form_lhs = "outcome")

  # Perform assertions on the result
  expect_message(smdi_outcome(data = data, model = "glm", form_lhs = "outcome"))
  expect_true("covariate4" %in% result$covariate)
  expect_equal(nrow(result), 1)
})

# Test case for linear model with continuous outcome
test_that("Linear model with continuous outcome works correctly", {
  set.seed(42)
  data <- data.frame(outcome = rnorm(100),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  result <- smdi_outcome(data = data, model = "linear", form_lhs = "outcome")

  # Perform assertions on the result
  expect_true("covariate4" %in% result$covariate)
  expect_equal(nrow(result), 1)
})

# Test case for Cox model with time-to-event outcome
test_that("Cox model with time-to-event outcome works correctly", {
  library(survival)
  set.seed(42)
  data <- data.frame(time =  rexp(100, rate = 0.2),
                     event = rbinom(100, 1, 0.5),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  result <- smdi_outcome(data = data, model = "cox", form_lhs = "Surv(time, event)")

  # Perform assertions on the result
  expect_true("covariate4" %in% result$covariate)
  expect_equal(nrow(result), 1)
})

# error with glm family
test_that("<glm_family> is specified although <model> is either 'linear' or 'cox'", {
  set.seed(42)
  data <- data.frame(outcome = rbinom(100, 1, 0.5),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  # Perform assertions on the result
  expect_error(smdi_outcome(data = data, model = "linear", glm_family = binomial(link = 'logit'), form_lhs = "outcome"))
})

# missings in outcome variable
test_that("NAs in outcome variable", {
  data <- data.frame(outcome = c(NA, NA, rbinom(98, 1, 0.5)),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  # Perform assertions on the result
  expect_error(smdi_outcome(data = data, model = "glm", glm_family = binomial(link = 'logit'), form_lhs = "outcome"))
})

# lifecycle error
test_that("model='logistic' lifecycle warning", {
  set.seed(42)
  data <- data.frame(outcome = rbinom(100, 1, 0.5),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  # Perform assertions on the result
  expect_error(smdi_outcome(data = data, model = "logistic", form_lhs = "outcome"))
})

# glm family
test_that("Incompatible model and glm_family calls", {
  set.seed(42)
  data <- data.frame(outcome = rbinom(100, 1, 0.5),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  # Perform assertions on the result
  expect_error(smdi_outcome(data = data, model = "cox", glm_family =  poisson(link = "log"), form_lhs = "outcome"))
})

# glm family
test_that("No error for compatible model and glm_family calls", {
  set.seed(42)
  data <- data.frame(outcome = rbinom(100, 1, 0.5),
                     covariate1 = rnorm(100),
                     covariate2 = rnorm(100),
                     covariate3 = rnorm(100),
                     covariate4 = sample(c(1, 2, NA), 100, replace = TRUE))

  # Perform assertions on the result
  result_smdi <- smdi_outcome(data = data, model = "glm", glm_family =  poisson(link = "log"), exponentiated = F, form_lhs = "outcome")
  result_manual <- glm(formula = outcome ~ covariate1 + covariate2 + covariate3 + is.na(covariate4), data = data, family = poisson(link = "log"))

  expect_no_error(result_smdi)
  expect_contains(class(result_smdi), "data.frame")

  })
