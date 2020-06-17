# V3M-prediction

### R scripts
* `preprocessing_summary.R`: pre-processing and descriptive analyses. The resulting dataset and objects are contained in `data_incl_censor.RData`.
* `50fold_MI.R`: conducts 50-fold multiple imputation. The imputed datasets (and other objects) are stored in `imputed_50fold_assess6_incl_censor.RData`.
* `separate_aki_nonaki.R`: random forest and lasso models separating AKI and non-AKI subjects. Also pool the results from each imputed dataset.

### csv and RData files (available on OneDrive)
* `assess_composite_v6.csv`: this is where I store the composite survival *outcome*. The *covariates* in this file are not up-to-date and should *not* be used for analysis.
* `data_incl_censor.RData` and `imputed_50fold_assess6_incl_censor.RData`: as described above

For the purpose of re-doing the analysis (possibly with some modification on the RF/lasso models), directly loading `imputed_50fold_assess6_incl_censor.RData` and running `separate_aki_nonaki.R` should suffice.
