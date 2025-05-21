get_hyperpars <- function(m, p)
{
  
  del_mu = 1; mr_mu = 0
  del_sig2 = rep(1,p); mr_sig2 = rep(0,p)
  gel_alpha = rep(1,m); mr_alpha = rep(0,m) 
  W = W <- lhs::randomLHS(n = m, k = p)  
  
  list(del_mu = del_mu, mr_mu = mr_mu, del_sig2 = del_sig2, mr_sig2 = mr_sig2,
       del2_U = del2_U,  del_alpha = del_alpha, mr_alpha = mr_alpha, W = W)
}

get_initial_values <- function(m, p)
{
  mu = 0
  sig2 = rep(1,p)
  alpha = rep(1,m)
  list(mu = mu, alpha = alpha, sig2 = sig2)
}

YPPP <- function(Yna, na_indices, K, n_sams, n_burn, n_skip, m, p, n)
{
  
  I <- get_I(Yna) 
  ## hyper-parameters
  hyps      <- get_hyperpars(m, p)
  del_mu    <- hyps$del_mu
  mr_mu     <- hyps$mr_mu
  del_sig2  <- hyps$del_sig2
  mr_sig2   <- hyps$mr_sig2
  del2_U    <- hyps$del2_U
  del_alpha <- hyps$del_alpha
  mr_alpha  <- hyps$mr_alpha
  W         <- hyps$W
  
  ## initial values
  tmp   <- get_initial_values(m, p)
  mu    <- tmp$mu
  sig2  <- tmp$sig2
  alpha <- tmp$alpha
  
  ## posterior predictive probabilities
  y_ppp <- rep(0, sum(na_indices))
  ## chains
  B <- n_burn + n_skip*n_sams
  for (b in 1:B) 
  {
    Yna        <- sample_Y(I, na_indices, Yna, X, W, p, m, alpha, sig2, mu)
    tmp        <- sample_mu(b, del_mu, mr_mu, n, p, m, mu, alpha, sig2, W, X, Yna)
    mu         <- tmp$mu
    del_mu     <- tmp$del_mu
    mr_mu      <- tmp$mr_mu
    tmp        <- sample_sig2(b, del_sig2, mr_sig2, n, p, m, mu, alpha, sig2, W, X, Yna)
    sig2       <- tmp$sig2
    del_sig2   <- tmp$del_sig2
    mr_sig2    <- tmp$mr_sig2
    tmp        <- sample_alpha(b, del_alpha, mr_alpha, n, p, m, mu, alpha, sig2, W, X, Yna)
    alpha      <- tmp$alpha
    del_alpha  <- tmp$del_alpha
    mr_alpha   <- tmp$mr_alpha
    
    # ppp
    if ((b > n_burn) & (b%%n_skip == 0)) y_ppp <- y_ppp + Yna[na_indices]/n_sams
  }
  y_ppp
}


