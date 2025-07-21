rm(list = ls(all.names = TRUE))

loc <- 'raftery/'
base::source(paste0(loc, 'utils.R'))
Rcpp::sourceCpp(paste0(loc, 'cpp/cfunctions.cpp'))

seed <- 1
set.seed(seed)
p <- 2
n <- 40
m <- 10 * p
mu <- -.5
X <- matrix(runif(n * p), n, p)
g <- function(x, y) 1.5 * exp(x^2) * x^2

dat <- data_gen2(n, mu, X, g, T)
z <- dat$z # scores
Theta <- dat$Theta # probabilities
Yadj <- dat$Y # adjacency matrix
Y <- as.matrix(Yadj[lower.tri(Yadj)]) # vec upper triangular adjacency matrix

# graph
g <- igraph::graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F)
# density
den <- igraph::graph.density(g)
# degree distribution
deg <- igraph::degree(g)

# control-to-case rate
# RNHY pick r in a way that n0 = 50
md <- mean(deg)
nc <- 25
r <- round(nc / md)
n0 <- round(r * md)
n0 <- min(20, n0)

# strata
# undirected network: upper triangular matrix
# max stratum_labels: Dij = Inf
stratum_labels <- vector("list", n)
stratum_sizes <- vector("list", n)
stratum_names <- vector("list", n)
for (i in 1:(n - 1)) {
   sp <- c(igraph::shortest.paths(graph = g, v = i))
   sp <- sp[-c(1:i)]
   inf_indices <- sp == "Inf"
   sp[inf_indices] <- 1
   sp[inf_indices] <- max(sp) + 1
   stratum_labels[[i]] <- sp
   stratum_sizes[[i]] <- table(sp)
   stratum_names[[i]] <- as.numeric(names(stratum_sizes[[i]]))
   rm(sp, inf_indices)
}

# sampling for pilot MCMC
# if shortest distance is 1, then dyad inlusion probability is 1
# at least one dyad in each stratum, for shortest distance >= 2
set.seed(seed)
samples_pilot <- vector("list", n)
for (i in 1:(n - 1)) {
   # forced inclusion y_ij = 1
   si <- ((i + 1):n)[which(stratum_labels[[i]] == 1)]
   # forced inclusion y_ij = 0
   ni <- stratum_names[[i]]
   ni <- ni[ni > 1]
   for (h in ni) {
      tmp <- ((i + 1):n)[which(stratum_labels[[i]] == h)]
      si <- c(si, tmp[sample.int(length(tmp), 1)])
   }
   # SRS for remaining actors s.t. y_ij = 0
   nn <- n0 + sum(stratum_labels[[i]] == 1) - length(si)
   si <- c(si, sample(x = (1:n)[-c(i, si)], size = nn, replace = F))
   samples_pilot[[i]] <- sort(si)
   rm(tmp, ni, si, nn)
}
# accommodate samples upper triangular matrix
for (i in 1:(n - 1)) {
   for (j in (i + 1):n) {
      if (i %in% samples_pilot[[j]]) {
         samples_pilot[[i]] <- c(samples_pilot[[i]], j)
         samples_pilot[[j]] <- samples_pilot[[j]][which(
            samples_pilot[[j]] != i
         )]
      }
   }
   samples_pilot[[i]] <- sort(unique(samples_pilot[[i]]))
}

# expansion factors and strata
ns1 <- NULL
ns0 <- NULL
sa1 <- vector("list", n) # 1-based indices indices
sa0 <- vector("list", n) # 1-based indices indices
ef0 <- vector("list", n)
st0 <- vector("list", n)
for (i in 1:n) {
   if (length(samples_pilot[[i]]) > 0) {
      sai <- samples_pilot[[i]]
      lai <- stratum_labels[[i]][((i + 1):n) %in% sai]
      sai1 <- sai[lai == 1]
      sai0 <- sai[lai > 1]
      lai0 <- lai[lai > 1]
      ns1[i] <- length(sai1)
      ns0[i] <- length(sai0)
      Ni <- stratum_sizes[[i]][stratum_names[[i]] > 1]
      ni <- table(factor(x = lai0, levels = names(Ni)))
      nni0 <- as.numeric(names(ni))
      efi0 <- rep(NA, ns0[i])
      sti0 <- rep(NA, ns0[i])
      for (h in nni0) {
         efi0[lai0 == h] <- (Ni / ni)[nni0 == h]
         sti0[lai0 == h] <- h
      }
      sa1[[i]] <- as.numeric(sai1)
      sa0[[i]] <- as.numeric(sai0)
      ef0[[i]] <- as.numeric(efi0)
      st0[[i]] <- as.numeric(sti0)
      rm(sai, lai, sai1, sai0, lai0, efi0, sti0, nni0)
   } else {
      ns1[i] <- 0
      ns0[i] <- 0
      sa1[[i]] <- numeric(0)
      sa0[[i]] <- numeric(0)
      ef0[[i]] <- numeric(0)
      st0[[i]] <- numeric(0)
   }
}

# vectorize
sai1 <- NULL # 1-based indices indices
sai0 <- NULL # 1-based indices indices
efi0 <- NULL
for (i in 1:n) {
   sai1 <- c(sai1, sa1[[i]])
   sai0 <- c(sai0, sa0[[i]])
   efi0 <- c(efi0, ef0[[i]])
}

# pilot MCMC
model <- "LSSPh_pilot"
dataset <- "synthetic"
prefix <- paste0(model, "_", dataset)
path_outs <- loc
n_sams <- 10000
n_burn <- 50000
n_skip <- 1
set.seed(seed)
W <- lhs::randomLHS(n = m, k = p)
#MCMC_hat(sai1-1, sai0-1, ns1, ns0, efi0, Y, X, W, n, p, m, n_sams, n_burn, n_skip, prefix, path_outs)
out_mu <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_mu_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   1,
   T
)
out_sig2 <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_sig2_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   p,
   T
)
out_alpha <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_alpha_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   m,
   T
)

# determine sample sizes using pilot MCMC
sample_sizes <- vector("list", n)
for (i in 1:(n - 1)) {
   st <- stratum_names[[i]]
   st <- st[st > 1]
   H <- length(st)
   wi <- rep(0, H)
   for (b in 1:n_sams) {
      tmp <- NULL
      for (h in 1:H) {
         tmp[h] <- abs(loglik_hat_ih_R(
            i,
            sai0 = sa0[[i]][st0[[i]] == st[h]],
            efi0 = ef0[[i]][st0[[i]] == st[h]],
            p,
            m,
            mu = out_mu[b, ],
            alpha = out_alpha[b, ],
            sig2 = out_sig2[b, ],
            W,
            X
         ))
      }
      if (sum(tmp) > 0) {
         wi <- wi + abs(tmp / sum(tmp)) / n_sams
      } else {
         wi <- wi + 0
      }
   }
   if (sum(wi) > 0) {
      sample_sizes[[i]] <- round(ns0[i] * abs(wi / sum(wi)))
      ss <- stratum_sizes[[i]]
      ss <- ss[as.numeric(names(ss)) > 1]
      diff <- ss - sample_sizes[[i]]
      for (h in 1:H) {
         if (diff[h] < 0) {
            sample_sizes[[i]][h] <- ss[h]
         }
      }
      rm(ss, diff)
   } else {
      sample_sizes[[i]] <- rep(0, H)
   }
   if (i %% floor(0.1 * n) == 0)
      cat("Working on actor ", i, " out of ", n, "\n", sep = "")
   rm(st, H, wi, tmp)
}

# sampling for full MCMC
set.seed(seed)
samples_full <- vector("list", n)
for (i in 1:(n - 1)) {
   # forced inclusion y_ij = 1
   si <- ((i + 1):n)[which(stratum_labels[[i]] == 1)]
   # y_ij = 0
   ni <- stratum_names[[i]]
   ni <- ni[ni > 1]
   for (j in 1:length(ni)) {
      h <- ni[j]
      tmp <- ((i + 1):n)[which(stratum_labels[[i]] == h)]
      si <- c(
         si,
         tmp[sample.int(
            n = length(tmp),
            size = sample_sizes[[i]][j],
            replace = F
         )]
      )
   }
   samples_full[[i]] <- sort(si)
   rm(tmp, ni, si)
}

# accommodate samples upper triangular matrix
for (i in 1:(n - 1)) {
   for (j in (i + 1):n) {
      if (i %in% samples_full[[j]]) {
         samples_full[[i]] <- c(samples_full[[i]], j)
         samples_full[[j]] <- samples_full[[j]][which(samples_full[[j]] != i)]
      }
   }
   samples_full[[i]] <- sort(unique(samples_full[[i]]))
}

# expansion factors and strata
ns1 <- NULL
ns0 <- NULL
sa1 <- vector("list", n) # 1-based indices indices
sa0 <- vector("list", n) # 1-based indices indices
ef0 <- vector("list", n)
st0 <- vector("list", n)
for (i in 1:n) {
   if (length(samples_full[[i]]) > 0) {
      sai <- samples_full[[i]]
      lai <- stratum_labels[[i]][((i + 1):n) %in% sai]
      sai1 <- sai[lai == 1]
      sai0 <- sai[lai > 1]
      lai0 <- lai[lai > 1]
      ns1[i] <- length(sai1)
      ns0[i] <- length(sai0)
      Ni <- stratum_sizes[[i]][stratum_names[[i]] > 1]
      ni <- table(factor(x = lai0, levels = names(Ni)))
      nni0 <- as.numeric(names(ni))
      efi0 <- rep(NA, ns0[i])
      sti0 <- rep(NA, ns0[i])
      for (h in nni0) {
         efi0[lai0 == h] <- (Ni / ni)[nni0 == h]
         sti0[lai0 == h] <- h
      }
      sa1[[i]] <- as.numeric(sai1)
      sa0[[i]] <- as.numeric(sai0)
      ef0[[i]] <- as.numeric(efi0)
      st0[[i]] <- as.numeric(sti0)
      rm(sai, lai, sai1, sai0, lai0, efi0, sti0, nni0)
   } else {
      ns1[i] <- 0
      ns0[i] <- 0
      sa1[[i]] <- numeric(0)
      sa0[[i]] <- numeric(0)
      ef0[[i]] <- numeric(0)
      st0[[i]] <- numeric(0)
   }
}

# vectorize
sai1 <- NULL # 1-based indices indices
sai0 <- NULL # 1-based indices indices
efi0 <- NULL
for (i in 1:n) {
   sai1 <- c(sai1, sa1[[i]])
   sai0 <- c(sai0, sa0[[i]])
   efi0 <- c(efi0, ef0[[i]])
}

# model fitting
loc <- 'C:/Users/David.solano/OneDrive - Ipsos/David/Trabjo de Grado Maestria/CODIGOS/Sampling Version/Synthetic_1/'
model <- "LSSPh_full"
dataset <- "synthetic1"
prefix <- paste0(model, "_", dataset)
path_outs <- loc
n_sams <- 10000
n_burn <- 50000
n_skip <- 10
set.seed(seed)
W <- lhs::randomLHS(n = m, k = p)
ptm <- proc.time()
#MCMC_hat(sai1-1, sai0-1, ns1, ns0, efi0, Y, X, W, n, p, m, n_sams, n_burn, n_skip, prefix, path_outs)
ptm <- proc.time() - ptm

### load samples
out_ll <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_ll_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   1,
   T
)
out_mu <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_mu_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   1,
   T
)
out_sig2 <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_sig2_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   p,
   T
)
out_alpha <- matrix(
   c(as.matrix(read.table(
      paste0(loc, prefix, "_alpha_chain.txt"),
      quote = "\"",
      comment.char = ""
   ))),
   n_sams,
   m,
   T
)

### effective sample sizes
summary(coda::effectiveSize(cbind(out_mu, out_sig2, out_alpha)))

chain_plot_plotly = function(x, main = '', truth = NULL) {
   data <- data.frame(x = 1:length(x), y = x)
   if (!is.null(truth))
      fig <- plot_ly(
         data,
         x = ~x,
         y = ~y,
         type = 'scatter',
         mode = 'lines',
         mode = 'lines'
      ) %>%
         add_lines(y = truth, line = list(color = 'rgb(205, 12, 24)')) %>%
         layout(
            title = main,
            font = 14,
            showlegend = FALSE,
            xaxis = list(showgrid = F),
            yaxis = list(showgrid = F)
         )

   if (is.null(truth))
      fig <- plot_ly(
         data,
         x = ~x,
         y = ~y,
         type = 'scatter',
         mode = 'lines',
         mode = 'lines'
      ) %>%
         layout(
            title = main,
            font = 14,
            showlegend = FALSE,
            xaxis = list(showgrid = F),
            yaxis = list(showgrid = F)
         )

   return(fig)
}

### convergence (chains)

chain_plot_plotly(x = out_ll, main = "Log-likelihood")
chain_plot_plotly(x = out_mu, main = 'mu', truth = mu)
for (i in 1:p)
   print(chain_plot_plotly(x = out_sig2[, i], main = paste0('sigma_', i)))
for (i in sample(1:m, 2))
   print(chain_plot_plotly(x = out_alpha[, i], main = paste0('alpha_', i)))


chain_plot(x = out_ll, main = "Log-likelihood")
chain_plot(x = out_mu, main = expression(mu), truth = mu)
hist_plot(x = out_ll, main = "Log-likelihood")
hist_plot(x = out_mu, main = expression(mu), truth = mu)
for (i in 1:p) hist_plot(x = out_sig2[, i], main = bquote(sigma[.(i)]^2))
for (i in 1:p) chain_plot(x = out_sig2[, i], main = bquote(sigma[.(i)]^2))
for (i in sample(1:m, 2))
   chain_plot(x = out_alpha[, i], main = bquote(alpha[.(i)]))

### surface
out_z <- matrix(NA, n_sams, n)
for (i in 1:n) {
   for (b in 1:n_sams) {
      out_z[b, i] <- zeta(p, m, X[i, ], out_alpha[b, ], out_sig2[b, ], W)
   }
}

### surface translation
out_z <- out_z - min(colMeans(out_z))

# posterior
for (i in sample(1:n, 2))
   hist_plot(x = out_z[, i], main = paste0("z_", i), truth = z[i])

figures = NULL

vline <- function(x = 0, color = "red") {
   list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color)
   )
}

for (i in sample(1:n, 2))
   figures[[paste0('z_', i)]] = plot_ly(
      x = out_z[, i],
      type = "histogram",
      alpha = 0.6,
      histnorm = 'probability',
      marker = list(color = 'rgb(158,202,225)', width = 1.5)
   ) %>%
      #add_lines(x = z[i], y = c(0, 0.05), line = list(color = 'rgb(205, 12, 24)')) %>%
      layout(
         shapes = list(vline(z[i])),
         title = "",
         xaxis = list(title = paste0('z_', i)),
         yaxis = list(title = ""),
         font = list(size = 14),
         showlegend = FALSE
      )
figures

# differences between true scores and estimated scores
round(summary(z - colMeans(out_z)), 3)

# estimated surface
x <- seq(from = 0, to = 1, length = m)
y <- seq(from = 0, to = 1, length = m)
z_hat <- matrix(0, m, m)
for (i in 1:m) {
   for (j in 1:m) {
      for (b in 1:n_sams) {
         z_hat[i, j] <- z_hat[i, j] +
            zeta(p, m, c(x[i], y[j]), out_alpha[b, ], out_sig2[b, ], W) / n_sams
      }
   }
}
nrz <- nrow(z_hat)
ncz <- ncol(z_hat)
nbcol <- 100
jet.colors <- colorRampPalette(c("blue", "green"))
color <- jet.colors(nbcol)
zfacet <- (z_hat[-1, -1] +
   z_hat[-1, -ncz] +
   z_hat[-nrz, -1] +
   z_hat[-nrz, -ncz]) /
   4
windows(width = 7.5, height = 5)
par(
   mfrow = c(2, 3),
   mar = c(3, 3, 2, 1),
   mgp = c(1.6, 0.6, 0),
   oma = c(0, 0, 0, 0)
)
for (theta in c(30, 60, 120, 150, 210, 240))
   persp(
      x = x,
      y = y,
      z = z_hat,
      theta = theta,
      phi = 30,
      col = color[cut(zfacet, nbcol)],
      axes = T,
      ticktype = "simple",
      cex.axis = 0.5,
      zlab = "z"
   )
title(main = "Estimated surface", outer = T, line = -1)

### interaction probabilities
Theta_hat <- matrix(0, n, n)
for (i in 1:(n - 1)) {
   for (j in (i + 1):n) {
      for (b in 1:n_sams) {
         Theta_hat[i, j] <- Theta_hat[i, j] +
            exp(pnorm(out_mu[b] - abs(out_z[b, i] - out_z[b, j]), log.p = T)) /
               n_sams
      }
   }
}
Theta_hat <- Theta_hat + t(Theta_hat)

# differences between true probabilities and estimated probabilities
round(summary(c(Theta) - c(Theta_hat)), 3)

windows(width = 7, height = 3.5)
par(
   mfrow = c(1, 2),
   mar = c(3, 3, 2, 1),
   mgp = c(1.6, 0.6, 0),
   oma = c(0, 0, 0, 0)
)
heat.plot0(Theta, main = "True probabilities")
heat.plot0(Theta_hat, main = "Estimated probabilities")

# 15984.89 secs = 4.44 hrs (150000 iterations)
# 0.1065 secs per iteration
ptm / (n_burn + n_sams * n_skip)

# Effective sample sizes
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3992    4478    5752    6010    7209    8596

# LSSP scores differences
# Estimates are underestimating true scores
# It seems that alphas are not fully identifiable due to rotations and
# reflections, but surface expansions and absolute differences are (see proof
# using orthogonal matrices).
# Surfaces are not identifiable due to translations. A partial solution consists
# in translating the surface's minimum to zero.

# Probabilities differences
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.157  -0.022  -0.006  -0.010   0.006   0.115

# in-sample model checking
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

#-------------------------------------------------------------------------------

test_stats <- c("Density", "Transitivity", "Assortativity")
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
   rm(Yrep, g)
   if (b %% (0.1 * n_sams) == 0)
      cat("Completed ", b / n_sams * 100, "% \n", sep = "")
}

# observed stats
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

for (i in 1:n_stats)
   hist_plot(
      x = out_test_stats[, i],
      main = test_stats[i],
      truth = obs_stats[i]
   )
