# smdi 0.3.0

* Incorporated comments from peer-review in JAMIA Open (Weberpals et al. 2024, <doi:10.1093/jamiaopen/ooae008>)

* Added `tune` parameter to `smdi_rf` to allow users to perform 5-fold cross validation and optimized random search for `mtry`()

* Changes to `smdi_outcome`: the `model` parameter option does not accept `logistic` anymore for logistic regressions but `glm` along with a new corresponding parameter `glm_family` to allow users to take advantage of all possible `glm` families as an outcome regression model (CAVE: no backwards compatibility)

* Variables are now one-hot encoded before running *naniar::mcar_test()* in `smdi_little` to address potential issues with categorical variables and to be consistent with smdi_hotelling. Results may differ slightly from those form previous versions and we suggest re-running analyses.

* Changed n_cores from a warning to a message notifying the user

* Improvement to `smdi_style_gt` to show correct formatting in `gt` exports of any supported type

* General maintenance and dependency management

# smdi 0.2.2

* CRAN release 

* Formally implemented unit tests

* Added unit test coverage report to pkgdown website

* Implemented automated GitLab CI/CD pipeline to run checks on daily basis

* Minor fixes and improvements in documentation of functions


# smdi 0.2.1

* Included re-exports of naniar's `gg_miss_upset` and mice's `md.pattern` functions to explore missing data patterns.

* New function `smdi_style_gt()` to make publication-ready tables based on objects of class smdi in combination with the `gt()` package.

* Added more details to *Routine structural missing data diagnostics* vignette.

* Updated `README` with more details and guidance on how to interpret the three group diagnostics and apply those to a real-world study.

* Some improved documentation here and there.

# smdi 0.2.0

* `smdi_asmd()`, and consequently also `smdi_diagnose()`, now also outputs the minimum (min) and maximum (max) absolute standardized mean difference (asmd) in addition to the mean/median to provide more comprehensive information about the asmd range without having to look at each asmd plot individually.

* In case of monotone missing data patterns, we observed unreasonably high AUC values for the Group 2 diagnostic which was caused by other partially observed covariates being almost perfect linear predictors of missingness. The new version has an in-built mechanism to prompt a message if AUCs are very high (> 0.9). The prompt also gives additional details about the covariate for which this behavior was observed and the strongest predictor based on the mean decrease in accuracy. In case of monotonicity this typically another partially observed covariate which would then be flagged with a “_NA” suffix. Based on the prompt, the analyst can then decide if this variable should be better dropped for the smdi diagnostics.

* To address issues and learning around multivariate missing data and handling of monotone missing data patterns in `smdi`, an additional vignette on `Multivariate missingness and monotonicity` was added.

* Change of colors in plots produced by `smdi_rf()` to address color-blindness

* Some improved documentation for smdi_diagnose, smdi_asmd and smdi_rf

# smdi 0.1.0

* Internal release of version 0.1.0 for beta testing

* First draft of all `smdi_xxx()` functions.

* Implementation of parallel processing to increase computational speed using `mclapply` (UNIX machines only)

* Initial build of website using `pkgdown`.

* Added a `NEWS.md` file to track changes to the package.

* Created three vignettes to learn more about the `smdi` package


