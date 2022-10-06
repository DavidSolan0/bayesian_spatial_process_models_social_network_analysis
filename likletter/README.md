This folder contains the codes implemented to do a complete analysis of the network data. Therefore, we strongly encourage the reader to use them in the following order:

* [descriptive](https://github.com/DavidSolan0/bayesian_spatial_process_models_social_network_analysis/blob/main/likletter/fit/descriptive.R): Here we do a data **descriptive** analysis. In case that you more than two covariates and want to study the estimated surface this EDA can give you insights to determine which couple of covariates you should conserv. 

* [fit](https://github.com/DavidSolan0/bayesian_spatial_process_models_social_network_analysis/blob/main/likletter/fit/fit.R): Here we charge the *Rcpp* codes with the MCMC algorithm and fit the model, keep in mind that this could be a long process. To train the model, we recommend rescaling the covariates to the $[0,1]$ interval. 

* [gof](https://github.com/DavidSolan0/bayesian_spatial_process_models_social_network_analysis/blob/main/likletter/fit/gof.R): Here we do a goodness-of-fit evaluation of the model with specific metrics of interest. You could use these metrics to test your model's performance or use more suitable ones if that is the case. 



