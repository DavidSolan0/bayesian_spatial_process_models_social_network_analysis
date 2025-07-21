loc <- "likletter/"
model <- "LSSP"
dataset <- "synthetic"
prefix <- paste0(model, "_", dataset)
path_outs <- paste0(loc, "results/fit/")
seed <- 1

n_sams <- 10000
p <- 2
n <- 40
m <- 10 * p

base::source(paste0(loc, "utils.R"))
Rcpp::sourceCpp(paste0(loc, "fit/src/cpp/cfunctions.cpp"))

### load samples
out_ll <- read_chain_file(path_outs, prefix, "ll", n_sams, 1)
out_mu <- read_chain_file(path_outs, prefix, "mu", n_sams, 1)
out_sig2 <- read_chain_file(path_outs, prefix, "sig2", n_sams, p)
out_alpha <- read_chain_file(path_outs, prefix, "alpha", n_sams, m)

### Data
set.seed(seed)
X <- matrix(runif(n * p), n, p)
g <- function(x, y) 1.5 * exp(x^2) * x^2
mu = -.5

dat <- data_gen2(n, mu, X, g)
z <- dat$z # scores
Theta <- dat$Theta # probabilities
Yadj <- dat$Y # adjacency matrix

set.seed(seed)
Y <- as.matrix(Yadj[lower.tri(Yadj)])
W <- lhs::randomLHS(n = m, k = p)

### surface
out_z <- matrix(NA, n_sams, n)
for (i in 1:n) {
  for (b in 1:n_sams) {
    out_z[b, i] <- zeta(p, m, X[i, ], out_alpha[b, ], out_sig2[b, ], W)
  }
}

### surface translation
out_z <- out_z - min(colMeans(out_z))

### in-sample model checking
test_stats <- c(
  "Density",
  "Transitivity",
  "Assortativity",
  'Diameter',
  'Avg_Distance',
  'Avg_Degree',
  'Sd_Degree'
)
n_stats <- length(test_stats)
out_test_stats <- matrix(NA, n_sams, n_stats)
colnames(out_test_stats) <- test_stats
set.seed(1234)
for (b in 1:n_sams) {
  Yrep <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      Yrep[i, j] <- Yrep[j, i] <- rbinom(
        1,
        1,
        exp(pnorm(out_mu[b] - abs(out_z[b, i] - out_z[b, j]), log.p = T))
      )
    }
  }
  g <- graph.adjacency(adjmatrix = Yrep, mode = "undirected", diag = F)
  out_test_stats[b, 1] <- igraph::graph.density(g)
  out_test_stats[b, 2] <- igraph::transitivity(g)
  out_test_stats[b, 3] <- igraph::assortativity.degree(g)
  out_test_stats[b, 4] <- igraph::diameter(g)
  out_test_stats[b, 5] <- igraph::mean_distance(g)
  out_test_stats[b, 6] <- mean(igraph::degree(g))
  out_test_stats[b, 7] <- sd(igraph::degree(g))
  rm(Yrep, g)
  if (b %% (0.1 * n_sams) == 0)
    cat("Completed ", b / n_sams * 100, "% \n", sep = "")
}

### observed stats
obs_stats <- NULL
obs_stats[1] <- igraph::graph.density(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
))
obs_stats[2] <- igraph::transitivity(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
))
obs_stats[3] <- igraph::assortativity.degree(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
))
obs_stats[4] <- igraph::diameter(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
))
obs_stats[5] <- igraph::mean_distance(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
))
obs_stats[6] <- mean(igraph::degree(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
)))
obs_stats[7] <- sd(igraph::degree(graph.adjacency(
  adjmatrix = Yadj,
  mode = "undirected",
  diag = F
)))

figures = list()
for (i in c(1:3, 5:7))
  figures[[test_stats[i]]] = plot_ly(
    x = out_test_stats[, i],
    type = "histogram",
    alpha = 0.6,
    marker = list(color = 'rgb(158,202,225)', width = 1.5)
  ) %>%
    add_trace(
      x = c(obs_stats[i], obs_stats[i]),
      y = c(0, 500),
      type = "scatter",
      mode = "markers+lines",
      line = list(color = 'rgb(205, 12, 24)')
    ) %>%
    layout(
      title = "",
      xaxis = list(title = test_stats[i]),
      yaxis = list(title = ""),
      font = list(size = 14),
      showlegend = FALSE
    )
#hist_plot(x = out_test_stats[,i], main = test_stats[i], truth = obs_stats[i])

figures
subplot(figures[-4], nrows = 3)
