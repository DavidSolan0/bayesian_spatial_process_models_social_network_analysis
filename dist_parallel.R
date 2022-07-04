cv_parallel <- function(Y, K_range, n_burn, n_sams, n_skip, model, dataset, loc)
{
  ## info
  cat(" * CV", "\n", sep = "")
  ## indices
  L <- 5
  set.seed(1234)
  folds <- get_folds(M = nrow(Y), 1, L)
  KL <- as.matrix(expand.grid(K_range, 1:L))
  ## set cores and cluster
  suppressMessages(suppressWarnings(require(doParallel)))
  cl <- makeCluster(min(detectCores(), nrow(KL)))
  registerDoParallel(cl)
  cv <- foreach (i = 1:nrow(KL), .inorder = F, .packages = c("Rcpp")) %dopar% {
    #source   (paste0(loc, "code/", "rfunctions.R"))
    #source   (paste0(loc, "code/", model, "_functions.R"))
    #sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
    
    base::source(paste0(loc, "dist_functions.R"))
    Rcpp::sourceCpp(paste0(loc, "dist_functions.cpp"))
    
    K <- KL[i,1]
    l <- KL[i,2]
    na_indices <- folds == l
    y_true <- as.factor(Y[na_indices])
    Yna <- Y
    Yna[na_indices] <- NA
    y_ppp <- YPPP(Yna, na_indices, K, n_sams, n_burn, n_skip)
    cv <- list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
    save(cv, file = paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_l_", l, "_", dataset, "_", "cv.RData"))
  }
  stopCluster(cl)
}
