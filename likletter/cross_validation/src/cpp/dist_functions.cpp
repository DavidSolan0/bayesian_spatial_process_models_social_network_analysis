#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;


uword get_k (const int& i, const int& ii, const double& I)
{
  return(0.5*(I*(I-1.0) - (I-i)*(I-i-1.0)) + ii - i - 1.0);
} 

double round_to_digits (double value, int digits)
{
  // otherwise it will return 'nan' due to the log10() of zero
  if (value == 0.0) return 0.0;
  double factor = pow(10.0, digits - ceil(log10(fabs(value))));
  return round(value*factor)/factor;   
}

char* mypaste0 (string path, string prefix, string name)
{
  // pastes path and name together
  stringstream strname;
  strname << path << prefix << "_" << name << ".txt";
  string fullname = strname.str();
  string::iterator p = fullname.begin();
  char* chr = &(*p);
  return( chr );
}

double tunning (const double& del, const double& mix_rate)
{ 
  // tunning paramter calibration
  double tmp = del;
  if (abs(mix_rate - 0.375) > 0.075) {
    int cont = 0;
    do {
      cont++;
      tmp = del + (0.1/double(cont))*(mix_rate - 0.375);
    } while ((tmp <= 0.0) && (cont <= 100));
    tmp = abs(tmp);
  } 
  return(tmp);
}

//[[Rcpp::export]]
inline double zeta (const int& p, const int& m, const rowvec& x, const rowvec& alpha, const rowvec& sig2, const mat& W) 
{
  // LSSP score
  // W : m x p matrix of locations (knots)
  int r, d;
  double tmp0 = -0.5*double(p)*std::log(2.0*datum::pi), tmp1 = 0.0, out = 0.0;
  for (r = 0; r < m; r++) {
    tmp1 = 0.0;
    for (d = 0; d < p; d++)
      tmp1 -= std::pow(x.at(d) - W.at(r,d), 2.0)/(2.0*sig2.at(d));
      out += alpha.at(r)*std::exp(tmp0 + tmp1);
  }
  return(out);
}

double loglik (const int& n, const int& p, const int& m, const double& mu, const rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
  // log-likelihood
  int i, j;
  double eta, out = 0.0;
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      eta = mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(j), alpha, sig2, W));
      out += R::pnorm(eta, 0, 1, Y.at(0.5*(n*(n-1)-(n-i)*(n-i-1))+j-i-1), 1);
    }
  }
  return(out);
}

// [[Rcpp::export]]
List sample_mu (const int& b, double& del_mu, double& mr_mu, const int& n, const int& p, const int& m, double& mu, const rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
  // sample mu (metropolis-hastings step)
  int i, j, y;
  double mr = 0.0, mu_p, logr, ados;  // v2_mu = 10.0;
  // proposal
  mu_p = R::rnorm(mu, del_mu);
  // acceptance probability
  logr = (std::pow(mu_p, 2.0) - std::pow(mu, 2.0))/(-2.0*10.0);
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      y = Y.at(0.5*(n*(n-1)-(n-i)*(n-i-1))+j-i-1);
      ados = abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(j), alpha, sig2, W));
      logr += R::pnorm(mu_p - ados, 0, 1, y, 1) - R::pnorm(mu - ados, 0, 1, y, 1);
    }
  }
  // set
  if (R::runif(0,1) < std::exp(logr)) {
    mr++;
    mu = mu_p;
  }
  mr_mu = ((double(b)-1.0)*mr_mu + mr)/double(b);
  if (b % 100 == 0)
    del_mu = tunning(del_mu, mr_mu);
  
  return List::create(Named("mu")     = mu,
                      Named("del_mu") = del_mu,
                      Named("mr_mu")  = mr_mu);
  
}

// [[Rcpp::export]]
List sample_sig2 (const int& b, rowvec& del_sig2, rowvec& mr_sig2, const int& n, const int& p, const int& m, const double& mu, const rowvec& alpha, rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
  // sample sig2 (metropolis-hastings step)
  int d, i, j, y;
  double a_sig = 3.0, b_sig = (a_sig - 1.0)*(1.0/(2.0*std::log(2.0))); // prior
  double mr, sig2_c, sig2_p, the_c, the_p, logr;
  for (d = 0; d < p; d++) {
    mr = 0.0;
    // proposal
    sig2_c = sig2.at(d);
    the_c  = std::log(sig2_c);
    the_p  = R::rnorm(the_c, del_sig2.at(d));
    sig2_p = std::exp(the_p);
    // acceptance probability
    logr = -(a_sig + 1.0)*(std::log(sig2_p) - std::log(sig2_c)) - b_sig*(1.0/sig2_p - 1.0/sig2_c); // prior
    for (i = 0; i < n-1; i++) {
      for (j = i+1; j < n; j++) {
        y = Y.at(0.5*(n*(n-1)-(n-i)*(n-i-1))+j-i-1);
        sig2.at(d) = sig2_p;
        logr += R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(j), alpha, sig2, W)), 0, 1, y, 1);
        sig2.at(d) = sig2_c;
        logr -= R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(j), alpha, sig2, W)), 0, 1, y, 1);
      }
    }
    logr += the_p - the_c; // jacobian
    // set
    if (R::runif(0,1) < std::exp(logr)) {
      mr++;
      sig2.at(d) = sig2_p;
    } else {
      sig2.at(d) = sig2_c;
    }
    mr_sig2.at(d) = ((double(b)-1.0)*mr_sig2.at(d) + mr)/double(b);
    if (b % 100 == 0)
      del_sig2.at(d) = tunning(del_sig2.at(d), mr_sig2.at(d));
  }
  return List::create(Named("sig2")    = sig2,
                     Named("del_sig2") = del_sig2,
                     Named("mr_sig2")  = mr_sig2);
  
  
}

// [[Rcpp::export]]
List  sample_alpha (const int& b, rowvec& del_alpha, rowvec& mr_alpha, const int& n, const int& p, const int& m, const double& mu, rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
  // sample alpha (metropolis-hastings step)
  int d, i, j, y;
  double mr, alpha_c, alpha_p, logr; // v2_mu = 10.0;
  for (d = 0; d < m; d++) {
    mr = 0.0;
    // proposal
    alpha_c = alpha.at(d);
    alpha_p = R::rnorm(alpha_c, del_alpha.at(d));
    // acceptance probability
    logr = (pow(alpha_p, 2.0) - pow(alpha_c, 2.0))/(-2.0*10.0); // prior
    for (i = 0; i < n-1; i++) {
      for (j = i+1; j < n; j++) {
        y = Y.at(0.5*(n*(n-1)-(n-i)*(n-i-1))+j-i-1);
        alpha.at(d) = alpha_p;
        logr += R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(j), alpha, sig2, W)), 0, 1, y, 1);
        alpha.at(d) = alpha_c;
        logr -= R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(j), alpha, sig2, W)), 0, 1, y, 1);
      }
    }
    // set
    if (R::runif(0,1) < std::exp(logr)) {
      mr++;
      alpha.at(d) = alpha_p;
    } else {
      alpha.at(d) = alpha_c;
    }
    mr_alpha.at(d) = ((double(b)-1.0)*mr_alpha.at(d) + mr)/double(b);
    if (b % 100 == 0)
      del_alpha.at(d) = tunning(del_alpha.at(d), mr_alpha.at(d));
  }
  return List::create(Named("alpha")     = alpha,
                      Named("del_alpha") = del_alpha,
                      Named("mr_alpha")  = mr_alpha);
}

// [[Rcpp::export]]
vec sample_Y (const double& I, const vec& na_indices, vec Yna,
              const mat& X, const mat& W, const int& p, const int& m,
              const rowvec& alpha, const rowvec& sig2, double& mu)
{
  // sample NA values in Y
  uword k;
  double eta;
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      k = get_k(i, ii, I);
      eta = mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(ii), alpha, sig2, W));
      if (na_indices[k] == true) Yna[k] = R::rbinom(1, R::pnorm(eta,0,1,1,0));
    }
  }
  return(Yna);
}