# Bayesian Spatial Process Models for Social Network Analysis

This directory contains the implementation of Bayesian spatial process models for social network analysis.

## Directory Structure

- `src/`: Source code
  - `cpp/`: C++ implementation of core model functions
  - `R/`: R implementation
    - `model/`: Core model functions
    - `analysis/`: Analysis and visualization functions
- `scripts/`: Execution scripts
  - `train.R`: Model training script
  - `evaluate.R`: Model evaluation script
- `tests/`: Model testing and validation
  - `gof.R`: Goodness of fit tests
- `data/`: Data storage (if needed)
- `results/`: Generated results
  - `figures/`: Generated plots and visualizations
  - `models/`: Saved model outputs

## Usage

1. Training the model:
```R
Rscript scripts/train.R
```

2. Evaluating the model:
```R
Rscript scripts/evaluate.R
```

3. Running goodness of fit tests:
```R
Rscript tests/gof.R
```

## Dependencies

- R packages: igraph, plotly, lhs
- C++ dependencies: Rcpp

## Model Description

[Add brief description of the model and its components]

This corresponds to the computational implementation of the algorithm in the B. appendix in this [documment](https://repositorio.unal.edu.co/handle/unal/82234). 

