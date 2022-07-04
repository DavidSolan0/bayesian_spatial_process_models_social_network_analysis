rm(list = ls())

loc       <- "C:/Users/David.solano/OneDrive - Ipsos/David/Trabjo de Grado Maestria/CODIGOS/Trabajo de Grado - Codigos/Synthetic1/"
model     <- "LSSP"
dataset   <- "synthetic"
prefix    <- paste0(model,"_",dataset)
path_outs <- loc
seed      <- 1     

base::source   (paste0(loc, "rfunctions.R"))
Rcpp::sourceCpp(paste0(loc, "cfunctions.cpp"))

### data generation
### p: no. covariates
### n: no. actors
### m: no. knots (locations)
### X: covariates matrix (n x p)
### g: surface

set.seed(seed)
p  <- 2    
n  <- 40  
m  <- 10*p
mu <- -.5
X  <- matrix(runif(n*p), n, p)
g  <- function(x, y) 1.5*exp(x^2)*x^2

dat   <- data_gen2(n, mu, X, g)
z     <- dat$z      # scores
Theta <- dat$Theta  # probabilities 
Yadj  <- dat$Y      # adjacency matrix 
#(Graph <- dat$graph) # Interactive Graph
#(surface <- dat$surface) # Interactive surface

### model fitting
set.seed(seed)
Y <- as.matrix(Yadj[lower.tri(Yadj)])  # vec upper triangular adjacency matrix
W <- lhs::randomLHS(n = m, k = p)      # knots matrix (m x p)
n_sams <- 10000
n_burn <- 50000
n_skip <- 10

ptm <- proc.time()
set.seed(1234)
#MCMC(Y, X, W, n, p, m, n_sams, n_burn, n_skip, prefix, path_outs)
ptm <- proc.time() - ptm

### load samples
out_ll    <- matrix(c(as.matrix(read.table(paste0(loc, prefix,    "_ll_chain.txt"), quote="\"", comment.char=""))), n_sams, 1, T)
out_mu    <- matrix(c(as.matrix(read.table(paste0(loc, prefix,    "_mu_chain.txt"), quote="\"", comment.char=""))), n_sams, 1, T)
out_sig2  <- matrix(c(as.matrix(read.table(paste0(loc, prefix,  "_sig2_chain.txt"), quote="\"", comment.char=""))), n_sams, p, T)
out_alpha <- matrix(c(as.matrix(read.table(paste0(loc, prefix, "_alpha_chain.txt"), quote="\"", comment.char=""))), n_sams, m, T)

### effective sample sizes
summary(coda::effectiveSize(cbind(out_mu, out_sig2, out_alpha)))

chain_plot_plotly = function(x, main = '', truth = NULL){
  
  data <- data.frame(x = 1:length(x), y = x)
  if(!is.null(truth)) 
    fig <- plot_ly(data, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                   mode = 'lines') %>% 
      add_lines(y = truth,line = list(color = 'rgb(205, 12, 24)')) %>% 
      layout(title = main, font = 14, showlegend = FALSE,
             xaxis = list(showgrid = F),
             yaxis = list(showgrid = F)) 
  
  if(is.null(truth)) 
    fig <- plot_ly(data, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                   mode = 'lines') %>% 
      layout(title = main, font = 14, showlegend = FALSE,
             xaxis = list(showgrid = F),
             yaxis = list(showgrid = F)) 
  
  return(fig)
}

### convergence (chains)
chain_plot_plotly(x = out_ll, main = "Log-likelihood")
chain_plot_plotly(x = out_mu, main = 'mu', truth = mu)
for (i in 1:p)
  print(chain_plot_plotly(x = out_sig2[,i],  main = paste0('sigma_',i)))
for (i in sample(1:m, 2))
  print(chain_plot_plotly(x = out_alpha[,i], main = paste0('alpha_',i)))

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
for (i in sample(1:n, 2))
  hist_plot(x = out_z[,i], main = paste0("z_",i), truth = z[i])

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

for(i in sample(1:n, 2))
  figures[[paste0('z_',i)]] = plot_ly(x = out_z[,i],
                                     type = "histogram",alpha = 0.6,histnorm='probability',
                                     marker = list(color = 'rgb(158,202,225)',width = 1.5)) %>%
  #add_lines(x = z[i], y = c(0, 0.05), line = list(color = 'rgb(205, 12, 24)')) %>% 
  layout(shapes = list(vline(z[i])), title = "",  xaxis = list(title = paste0('z_',i)),
         yaxis = list(title = ""),
         font = list(size = 14),showlegend = FALSE)
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

# differences between true probabilities and estimated probabilities
round(summary(c(Theta) - c(Theta_hat)), 3)


windows(width = 7, height = 3.5)
par(mfrow = c(1,2), mar = c(3,3,2,1), mgp = c(1.6,0.6,0), oma = c(0,0,0,0))
heat.plot0(Theta,     main = "True probabilities")
heat.plot0(Theta_hat, main = "Estimated probabilities")

# 15984.89 secs = 4.44 hrs (150000 iterations)
# 0.1065 secs per iteration
ptm/(n_burn + n_sams*n_skip)

# Effective sample sizes
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3587    3925    5184    5531    6950   10000 

# LSSP scores differences
# Estimates are underestimating true scores
# It seems that alphas are not fully identifiable due to rotations and 
# reflections, but surface expansions and absolute differences are (see proof 
# using orthogonal matrices).
# Surfaces are not identifiable due to translations. A partial solution consists
# in translating the surface's minimum to zero. 

# Probabilities differences
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.148  -0.019  -0.003  -0.006   0.011   0.141