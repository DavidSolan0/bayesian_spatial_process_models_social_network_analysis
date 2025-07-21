# raftery

This folder provides an alternative implementation of the Bayesian spatial process model for networks, following the approach of Raftery et al. (2012). It is similar in purpose to the main `likletter` folder, but with some key differences:

- **Computational Efficiency:**
  - The original Linkletter proposal (see `likletter/`) has computational cost $O(n^2)$ for a network of $n$ nodes, making it impractical for large networks.
  - Raftery et al. (2012) replace the full likelihood in the MCMC procedure with an unbiased estimation using an epidemiological case-control approach, reducing the cost to $O(n)$.

- **Implementation Style:**
  - The code here is less modular than in `likletter`, with more steps combined in single scripts.
  - Additional steps specific to the Raftery approach are included.

## Structure

- **src/R/model/fit.R**: Main script for fitting the model using the Raftery case-control likelihood approach.
- **src/cpp/cfunctions.cpp**: C++ code for computationally intensive routines, called from R.
- **utils.R**: Utility functions used throughout the analysis.
- **synthetic.R** (if present): Example or main script implementing the full workflow on synthetic data.

## How to Use

1. Edit and run the main R scripts (e.g., `src/R/model/fit.R` or `synthetic.R`) to fit the model and analyze results.
2. Adjust parameters and data paths as needed for your use case.
3. Results and outputs will be saved in the appropriate results directory.

## References

- Theoretical details can be found in Chapter 6 of this [document](https://repositorio.unal.edu.co/handle/unal/82234).
- For a more modular and general implementation, see the main [`likletter/`](../likletter/) folder.
