# cross_validation

This folder contains the code for performing cross-validation (CV) of the Bayesian spatial process model.

## Overview

- Implements a 5-fold cross-validation by default, but you can adjust the number of folds to suit your needs.
- Includes R and C++ scripts for efficient, parallelized computation of distance functions and model evaluation.

## Structure

- **src/R/**
  - `cv.R`: Main script to run cross-validation, including data partitioning, model training, and performance evaluation.
  - `dist_parallel.R`: Handles parallel computation of distance matrices. **Note:** Ensure the Rtools configuration matches your system.
  - `dist_functions.R`: Contains functions for distance calculations and imputation using the log-likelihood function.

- **src/cpp/**
  - `dist_functions.cpp`: C++ code for fast distance calculations, interfaced with R via Rcpp.

## How to Use

1. **Configure Parameters**  
   - By default, 5-fold CV is used. You can change the number of folds in `cv.R` as needed.

2. **Check Rtools**  
   - In `dist_parallel.R`, verify that the Rtools path and configuration are correct for your system.

3. **Run Cross-Validation**  
   - Execute `cv.R` to perform cross-validation and evaluate model predictive performance.

4. **Imputation and Log-Likelihood**  
   - For details on the imputation process using the log-likelihood function, refer to `dist_functions.R`.

## Output

- Results and performance metrics are saved in the `likletter/results/cv/` directory.

## References

- Theoretical details can be found in Chapter 2.2.4 of this [document](https://repositorio.unal.edu.co/handle/unal/82234).
- For a complete workflow, see the main [`likletter/README.md`](../README.md). 
