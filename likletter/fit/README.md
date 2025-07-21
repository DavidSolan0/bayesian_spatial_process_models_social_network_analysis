# fit

This folder contains the computational implementation of the Bayesian spatial process model algorithm, as described in Appendix B of this [document](https://repositorio.unal.edu.co/handle/unal/82234).

## Structure

- **src/R/model/**
  - `fit.R`: Main script for fitting the model using MCMC. Handles data preparation, model execution, and chain diagnostics.
  - `gof.R`: Computes goodness-of-fit metrics and performs posterior predictive checks.

- **src/R/analysis/**
  - `descriptive.R`: Performs exploratory and descriptive analysis of the network data, including degree distributions and network visualization.

- **src/cpp/**
  - `cfunctions.cpp`: C++ implementations of computationally intensive routines, interfaced with R via Rcpp.

## How to Use

1. **Descriptive Analysis**  
   Start with `src/R/analysis/descriptive.R` to explore your data and understand its structure.

2. **Model Fitting**  
   Use `src/R/model/fit.R` to fit the Bayesian spatial process model.  
   - Adjust parameters as needed for your dataset.
   - It is recommended to rescale covariates to the [0,1] interval for better model performance.

3. **Goodness-of-Fit**  
   Run `src/R/model/gof.R` to evaluate model fit using various network statistics and posterior predictive checks.

## Notes

- The scripts rely on utility functions from `likletter/utils.R`.
- C++ code is called via Rcpp for performance-critical tasks.
- Output files (e.g., MCMC chains) are saved in the `likletter/results/fit/` directory.

## References

- For theoretical background and algorithmic details, see Appendix B of this [document](https://repositorio.unal.edu.co/handle/unal/82234).
- For a complete workflow, refer to the main [`likletter/README.md`](../README.md). 

