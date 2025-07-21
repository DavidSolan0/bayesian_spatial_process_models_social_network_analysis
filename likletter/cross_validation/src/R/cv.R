rm(list = ls())

library(caret)
library(doParallel)
library(pROC)
library(ROCR)

loc <- "likletter/"
model <- "LSSP"
dataset <- "synthetic"
prefix <- paste0(model, "_", dataset)
path_outs <- paste0(loc, "results/cv/")
seed <- 1

get_folds <- function(M, J, L) {
  folds <- matrix(NA, M, J)
  for (j in 1:J) {
    for (l in 1:L) {
      indices <- which(is.na(folds[, j]))
      folds[sample(indices, floor(M / L), F), j] <- l
    }
    indices <- which(is.na(folds[, j]))
    if (sum(indices)) folds[indices, j] <- sample(1:L, 1, F)
  }
  folds
}

base::source(paste0(loc, "utils.R"))
base::source(paste0(loc, "cross_validation/src/R/dist_parallel.R"))

### data generation
### p: no. covariates
### n: no. actors
### m: no. knots (locations)
### X: covariates matrix (n x p)
### g: surface

set.seed(seed)
p <- 2
n <- 10
m <- 10 * p
mu <- -.5
X <- matrix(runif(n * p), n, p)
g <- function(x, y) 1.5 * exp(x^2) * x^2

dat <- data_gen2(n, mu, X, g)
z <- dat$z # scores
Theta <- dat$Theta # probabilities
Yadj <- dat$Y # adjacency matrix

### model fitting
Y <- as.matrix(Yadj[lower.tri(Yadj)]) # vec upper triangular adjacency matrix
W <- lhs::randomLHS(n = m, k = p) # knots matrix (m x p)
n_sams <- 10000
n_burn <- 50000
n_skip <- 10
K_range <- 2

cv_parallel(
  Y,
  X,
  K_range,
  m,
  p,
  n,
  n_burn,
  n_sams,
  n_skip,
  model,
  dataset,
  loc,
  path_outs
)

#* GOF
AUC = NULL
for (l in 1:5) {
  cv = read_(paste0(
    path_outs,
    "_K_",
    K_range,
    "_l_",
    l,
    "_",
    dataset,
    "_",
    "cv.RData"
  ))
  predictionsROCR = prediction(cv$y_ppp, cv$y_true)

  ROCRPerf = performance(predictionsROCR, "tpr", "fpr")
  plot(
    ROCRPerf,
    colorize = TRUE,
    print.cutoffs.at = seq(0, 1, 0.1),
    text.adj = c(-0.2, 1.7)
  )

  auc = c(auc, as.numeric(performance(predictionsROCR, "auc")@y.values))
  AUC = c(AUC, auc)
}

AUC = NULL
x11()
par(mfrow = c(3, 2))

for (l in 1:5) {
  load(paste0(
    path_outs,
    "_K_",
    K_range,
    "_l_",
    l,
    "_",
    dataset,
    "_",
    "cv.RData"
  ))

  predictionsROCR = prediction(cv$y_ppp, cv$y_true)

  ROCRPerf = performance(predictionsROCR, "tpr", "fpr")
  plot(
    ROCRPerf,
    colorize = TRUE,
    print.cutoffs.at = seq(0, 1, 0.1),
    text.adj = c(-0.2, 1.7)
  )

  auc = as.numeric(performance(predictionsROCR, "auc")@y.values)
  AUC = c(AUC, auc)
}

AUC
mean(AUC)
