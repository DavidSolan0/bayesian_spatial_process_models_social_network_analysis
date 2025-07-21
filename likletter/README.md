# likletter Folder Documentation

## Overview

The `likletter` folder contains all the code and utilities required to perform a complete Bayesian spatial process analysis on network data. The workflow is modular, with each subfolder and script handling a specific part of the analysis pipeline, from exploratory data analysis to model fitting, validation, and results visualization.

## Folder Structure

- **cross_validation/**
  - Implements cross-validation routines (default: 5-fold CV) for model assessment.
  - Contains R and C++ code for parallelized distance calculations and CV logic.
  - See `cross_validation/README.md` for details.

- **fit/**
  - Contains the main model fitting routines, including MCMC algorithms and goodness-of-fit checks.
  - R scripts for descriptive analysis, model fitting, and model evaluation.
  - C++ code for computationally intensive routines.
  - See `fit/README.md` for more.

- **results/**
  - Stores output files from model fitting and cross-validation, such as MCMC chains and performance metrics.

- **utils.R**
  - A collection of utility functions for plotting, data generation, and matrix operations used throughout the analysis.

## Key Scripts

- **fit/src/R/analysis/descriptive.R**
  - Performs exploratory and descriptive analysis of the network data, including degree distributions and network visualization.

- **fit/src/R/model/fit.R**
  - Main script for fitting the Bayesian spatial process model using MCMC. Handles data preparation, model execution, and chain diagnostics.

- **fit/src/R/model/gof.R**
  - Computes goodness-of-fit metrics and performs posterior predictive checks.

- **cross_validation/src/R/cv.R**
  - Orchestrates the cross-validation process, including data partitioning, model training, and performance evaluation.

- **cross_validation/src/R/dist_parallel.R** and **dist_functions.R**
  - Support scripts for distance calculations and parallel processing during cross-validation.

- **fit/src/cpp/cfunctions.cpp** and **cross_validation/src/cpp/dist_functions.cpp**
  - C++ implementations of computationally intensive functions, interfaced with R via Rcpp.

## How to Use

1. **Descriptive Analysis**:  
   Start with `fit/src/R/analysis/descriptive.R` to explore your data and understand its structure.

2. **Model Fitting**:  
   Use `fit/src/R/model/fit.R` to fit the Bayesian spatial process model. Adjust parameters as needed for your dataset.

3. **Goodness-of-Fit**:  
   Run `fit/src/R/model/gof.R` to evaluate model fit using various network statistics.

4. **Cross-Validation**:  
   Use `cross_validation/src/R/cv.R` to assess predictive performance via cross-validation.

5. **Results**:  
   Output files (e.g., MCMC chains) are saved in the `results/` directory for further analysis and visualization.

## References

- For theoretical background, see Chapter 2.2.4 and Appendix B of this [document](https://repositorio.unal.edu.co/handle/unal/82234).
- Each subfolder contains its own `README.md` with more specific instructions.

---

### Original Usage Order (Preserved)

This folder contains the codes implemented to do a complete analysis of the network data. Therefore, we strongly encourage the reader to use them in the following order:

* [descriptive](fit/src/R/analysis/descriptive.R): Here we do a data **descriptive** analysis. In case that you have more than two covariates and want to study the estimated surface, this EDA can give you insights to determine which couple of covariates you should conserve. 

* [fit](fit/src/R/model/fit.R): Here we load the *Rcpp* codes with the MCMC algorithm and fit the model. Keep in mind that this could be a long process. To train the model, we recommend rescaling the covariates to the $[0,1]$ interval. 

* [gof](fit/src/R/model/gof.R): Here we do a goodness-of-fit evaluation of the model with specific metrics of interest. You could use these metrics to test your model's performance or use more suitable ones if that is the case. 

* [cv](cross_validation/src/R/cv.R): Here we test the model predictive performance. 


