### Directorios y valores

loc       <- "C:/Users/David.solano/OneDrive - Ipsos/David/Trabjo de Grado Maestria/CODIGOS/Trabajo de Grado - Codigos/Synthetic1/"
model     <- "LSSP"
dataset   <- "synthetic"
prefix    <- paste0(model,"_",dataset)
path_outs <- loc
seed      <- 1  
n_sams <- 10000
p  <- 2    
n  <- 40  
m  <- 10*p

base::source   (paste0(loc, "rfunctions.R"))
Rcpp::sourceCpp(paste0(loc, "cfunctions.cpp"))

### load samples
out_ll    <- matrix(c(as.matrix(read.table(paste0(loc, prefix,    "_ll_chain.txt"), quote="\"", comment.char=""))), n_sams, 1, T)
out_mu    <- matrix(c(as.matrix(read.table(paste0(loc, prefix,    "_mu_chain.txt"), quote="\"", comment.char=""))), n_sams, 1, T)
out_sig2  <- matrix(c(as.matrix(read.table(paste0(loc, prefix,  "_sig2_chain.txt"), quote="\"", comment.char=""))), n_sams, p, T)
out_alpha <- matrix(c(as.matrix(read.table(paste0(loc, prefix, "_alpha_chain.txt"), quote="\"", comment.char=""))), n_sams, m, T)

### Data
set.seed(seed)
X  <- matrix(runif(n*p), n, p)
g  <- function(x, y) 1.5*exp(x^2)*x^2
mu = -.5

dat   <- data_gen2(n, mu, X, g)
z     <- dat$z      # scores
Theta <- dat$Theta  # probabilities 
Yadj  <- dat$Y      # adjacency matrix 

set.seed(seed)
Y <- as.matrix(Yadj[lower.tri(Yadj)]) 
W <- lhs::randomLHS(n = m, k = p) 

### effective sample sizes
summary(coda::effectiveSize(cbind(out_mu, out_sig2, out_alpha)))

### convergence (chains)
chain_plot(x = out_ll, main = "Log-likelihood")
chain_plot(x = out_mu, main = expression(mu), truth = mu)

for (i in 1:p)
  chain_plot(x = out_sig2[,i],  main = bquote(sigma[.(i)]^2))
for (i in sample(1:m, 2))
  chain_plot(x = out_alpha[,i], main = bquote(alpha[.(i)]))

### surface
out_z <- matrix(NA, n_sams, n)
for (i in 1:n) {
  for (b in 1:n_sams) {
    out_z[b,i] <- zeta(p, m, X[i,], out_alpha[b,], out_sig2[b,], W)
  }
}

### surface translation
out_z <- out_z - min(colMeans(out_z))

# posterior 
figures = list()
for (i in sample(1:n, 2))
  figures[[paste0("z_",i)]] = plot_ly(x = out_z[,i],
          type = "histogram",
          marker = list(color = 'rgb(158,202,225)', width = 1.5)) %>%
  add_trace(x = c(z[i],z[i]), y = c(0,500),
            type = "scatter", mode="markers+lines",
            line = list(color = 'rgb(205, 12, 24)')) %>% 
  layout(title = "",  xaxis = list(title = paste0("z_",i)),
         yaxis = list(title = ""),
         font = list(size = 14),showlegend = FALSE)

figures
#hist_plot(x = , main = paste0("z_",i), truth = z[i])

### estimated surface
x <- seq(from = 0, to = 1, length = m)
y <- seq(from = 0, to = 1, length = m)
z_hat <- matrix(0, m, m)
for (i in 1:m) {
  for (j in 1:m) {
    for (b in 1:n_sams) {
      z_hat[i,j] <- z_hat[i,j] + zeta(p, m, c(x[i],y[j]), out_alpha[b,], out_sig2[b,], W)/n_sams
    }
  }
} 
nrz <- nrow(z_hat)
ncz <- ncol(z_hat)
nbcol <- 100
jet.colors <- colorRampPalette(c("blue", "green"))
color <- jet.colors(nbcol)
zfacet <- (z_hat[-1, -1] + z_hat[-1, -ncz] + z_hat[-nrz, -1] + z_hat[-nrz, -ncz])/4
windows(width = 7.5, height = 5)
par(mfrow = c(2,3), mar = c(3,3,2,1), mgp = c(1.6,0.6,0), oma = c(0,0,0,0))
for (theta in c(30,60,120,150,210,240))
  persp(x = x, y = y, z = z_hat, theta = theta, phi = 30, col = color[cut(zfacet, nbcol)], 
        axes = T, ticktype = "simple", cex.axis = 0.5, zlab = "z")
title(main = "Estimated surface", outer = T, line = -1)

### interaction probabilities
Theta_hat <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    for (b in 1:n_sams) {
      Theta_hat[i,j] <- Theta_hat[i,j] + exp(pnorm(out_mu[b] - abs(out_z[b,i] - out_z[b,j]), log.p = T))/n_sams
    }
  }
}
Theta_hat <- Theta_hat + t(Theta_hat)

windows(width = 7, height = 3.5)
par(mfrow = c(1,1), mar = c(3,3,2,1), mgp = c(1.6,0.6,0), oma = c(0,0,0,0))
heat.plot0(Theta_hat, main = "Estimated probabilities")

### in-sample model checking
test_stats <- c("Density", "Transitivity", "Assortativity",'Diameter','Avg_Distance','Avg_Degree','Sd_Degree')
n_stats <- length(test_stats)
out_test_stats <- matrix(NA, n_sams, n_stats)
colnames(out_test_stats) <- test_stats
set.seed(1234)
for (b in 1:n_sams) {
  Yrep <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Yrep[i,j] <- Yrep[j,i] <- rbinom(1, 1, exp(pnorm(out_mu[b] - abs(out_z[b,i] - out_z[b,j]), log.p = T)))
    }
  }
  g <- graph.adjacency(adjmatrix = Yrep, mode = "undirected", diag = F)
  out_test_stats[b,1] <- igraph::graph.density(g)
  out_test_stats[b,2] <- igraph::transitivity(g)
  out_test_stats[b,3] <- igraph::assortativity.degree(g)
  out_test_stats[b,4] <- igraph::diameter(g)
  out_test_stats[b,5] <- igraph::mean_distance(g)
  out_test_stats[b,6] <- mean(igraph::degree(g))
  out_test_stats[b,7] <- sd(igraph::degree(g))
  rm(Yrep, g)
  if (b%%(0.1*n_sams)== 0) cat("Completed ", b/n_sams*100, "% \n", sep = "") 
}

### observed stats
obs_stats <- NULL
obs_stats[1] <- igraph::graph.density(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F))
obs_stats[2] <- igraph::transitivity(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F))
obs_stats[3] <- igraph::assortativity.degree(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F))
obs_stats[4] <- igraph::diameter(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F))
obs_stats[5] <- igraph::mean_distance(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F))
obs_stats[6] <- mean(igraph::degree(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F)))
obs_stats[7] <- sd(igraph::degree(graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F)))

figures = list()
for (i in c(1:3,5:7))
  figures[[test_stats[i]]] = plot_ly(x = out_test_stats[,i],
          type = "histogram",alpha = 0.6,
          marker = list(color = 'rgb(158,202,225)',width = 1.5)) %>%
  add_trace(x = c(obs_stats[i],obs_stats[i]), y = c(0,500),
            type = "scatter", mode="markers+lines",
            line = list(color = 'rgb(205, 12, 24)')) %>% 
  layout(title = "",  xaxis = list(title = test_stats[i]),
         yaxis = list(title = ""),
         font = list(size = 14),showlegend = FALSE)
#hist_plot(x = out_test_stats[,i], main = test_stats[i], truth = obs_stats[i])

figures
subplot(figures[-4],nrows=3)
