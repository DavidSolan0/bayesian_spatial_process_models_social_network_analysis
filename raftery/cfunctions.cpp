#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;
using namespace std;

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

//[[Rcpp::export]]
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

/*
//[[Rcpp::export]]
double loglik_hat_i (const int& i, const uvec& sai1, const uvec& sai0, const int& ni1, const int& ni0, const vec& efi0, const int& p, const int& m, const double& mu, const rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X)
{
        // log-likelihood hat i
        // sai1 : actors j s.t. y_ij = 1 (0-indexed based)
        // sai0 : actors j s.t. y_ij = 0 (0-indexed based)
        // ni1  : length(sai1)
        // ni0  : length(sai0)
        // efi0 : expansion factors of actors s.t. y_ij = 0
        double out = 0.0;
        for (int k = 0; k < ni1; k++) {
                out += R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[k]), alpha, sig2, W)), 0, 1, 1, 1);
        }
        for (int k = 0; k < ni0; k++) {
                out += efi0[k]*R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[k]), alpha, sig2, W)), 0, 1, 0, 1);
        }
        return(out);
}
*/

//[[Rcpp::export]]
double loglik_hat (const uvec& sai1, const uvec& sai0, const uvec& ns1, const uvec& ns0, const vec& efi0, const int& n, const int& p, const int& m, const double& mu, const rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X)
{
        // log-likelihood hat
        uword k;
        int i, i1 = 0, i0 = 0;
        double out = 0.0;
        for (i = 0; i < n; i++) {
                for (k = 0; k < ns1[i]; k++) {
                        out += R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[i1+k]), alpha, sig2, W)), 0, 1, 1, 1);
                }
                for (k = 0; k < ns0[i]; k++) {
                       out += efi0[i0+k]*R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[i0+k]), alpha, sig2, W)), 0, 1, 0, 1);
                }
                i1 += ns1[i];
                i0 += ns0[i];
        }
        return(out);
}

void sample_mu (const int& b, double& del_mu, double& mr_mu, const int& n, const int& p, const int& m, double& mu, const rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
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
}

void sample_mu_hat (const uvec& sai1, const uvec& sai0, const uvec& ns1, const uvec& ns0, const vec& efi0, const int& b, double& del_mu, double& mr_mu, const int& n, const int& p, const int& m, double& mu, const rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
        // sample mu (metropolis-hastings step)
        uword k;
        int i, i1 = 0, i0 = 0;
        double mr = 0.0, mu_p, logr, ados;  // v2_mu = 10.0;
        // proposal
        mu_p = R::rnorm(mu, del_mu);
        // acceptance probability
        logr = (std::pow(mu_p, 2.0) - std::pow(mu, 2.0))/(-2.0*10.0);
        for (i = 0; i < n; i++) {
                for (k = 0; k < ns1[i]; k++) {
                        ados  = abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[i1+k]), alpha, sig2, W));
                        logr += R::pnorm(mu_p - ados, 0, 1, 1, 1) - R::pnorm(mu - ados, 0, 1, 1, 1);
                }
                for (k = 0; k < ns0[i]; k++) {
                        ados  = abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[i0+k]), alpha, sig2, W));
                        logr += efi0[i0+k]*(R::pnorm(mu_p - ados, 0, 1, 0, 1) - R::pnorm(mu - ados, 0, 1, 0, 1));
                }
                i1 += ns1[i];
                i0 += ns0[i];
        }
        // set
        if (R::runif(0,1) < std::exp(logr)) {
                mr++;
                mu = mu_p;
        }
        mr_mu = ((double(b)-1.0)*mr_mu + mr)/double(b);
        if (b % 100 == 0)
                del_mu = tunning(del_mu, mr_mu);
}

void sample_sig2 (const int& b, rowvec& del_sig2, rowvec& mr_sig2, const int& n, const int& p, const int& m, const double& mu, const rowvec& alpha, rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
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
}

void sample_sig2_hat (const uvec& sai1, const uvec& sai0, const uvec& ns1, const uvec& ns0, const vec& efi0, const int& b, rowvec& del_sig2, rowvec& mr_sig2, const int& n, const int& p, const int& m, const double& mu, const rowvec& alpha, rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
        // sample sig2 (metropolis-hastings step)
        uword k;
        int d, i, i1, i0;
        double a_sig = 3.0, b_sig = (a_sig - 1.0)*(1.0/(2.0*std::log(2.0))); // prior
        double mr, sig2_c, sig2_p, the_c, the_p, logr;
        for (d = 0; d < p; d++) {
                i1 = 0, i0 = 0;
                mr = 0.0;
                // proposal
                sig2_c = sig2.at(d);
                the_c  = std::log(sig2_c);
                the_p  = R::rnorm(the_c, del_sig2.at(d));
                sig2_p = std::exp(the_p);
                // acceptance probability
                logr = -(a_sig + 1.0)*(std::log(sig2_p) - std::log(sig2_c)) - b_sig*(1.0/sig2_p - 1.0/sig2_c); // prior
                for (i = 0; i < n; i++) {
                        for (k = 0; k < ns1[i]; k++) {
                                sig2.at(d) = sig2_p;
                                logr += R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[i1+k]), alpha, sig2, W)), 0, 1, 1, 1);
                                sig2.at(d) = sig2_c;
                                logr -= R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[i1+k]), alpha, sig2, W)), 0, 1, 1, 1);
                        }
                        for (k = 0; k < ns0[i]; k++) {
                                sig2.at(d) = sig2_p;
                                logr += efi0[i0+k]*R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[i0+k]), alpha, sig2, W)), 0, 1, 0, 1);
                                sig2.at(d) = sig2_c;
                                logr -= efi0[i0+k]*R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[i0+k]), alpha, sig2, W)), 0, 1, 0, 1);
                        }
                        i1 += ns1[i];
                        i0 += ns0[i];
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
}

void sample_alpha (const int& b, rowvec& del_alpha, rowvec& mr_alpha, const int& n, const int& p, const int& m, const double& mu, rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
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
}

void sample_alpha_hat (const uvec& sai1, const uvec& sai0, const uvec& ns1, const uvec& ns0, const vec& efi0, const int& b, rowvec& del_alpha, rowvec& mr_alpha, const int& n, const int& p, const int& m, const double& mu, rowvec& alpha, const rowvec& sig2, const mat& W, const mat& X, const uvec& Y)
{
        // sample alpha (metropolis-hastings step)
        uword k;
        int d, i, i0, i1;
        double mr, alpha_c, alpha_p, logr; // v2_mu = 10.0;
        for (d = 0; d < m; d++) {
                i1 = 0, i0 = 0;
                mr = 0.0;
                // proposal
                alpha_c = alpha.at(d);
                alpha_p = R::rnorm(alpha_c, del_alpha.at(d));
                // acceptance probability
                logr = (pow(alpha_p, 2.0) - pow(alpha_c, 2.0))/(-2.0*10.0); // prior
                for (i = 0; i < n; i++) {
                        for (k = 0; k < ns1[i]; k++) {
                                alpha.at(d) = alpha_p;
                                logr += R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[i1+k]), alpha, sig2, W)), 0, 1, 1, 1);
                                alpha.at(d) = alpha_c;
                                logr -= R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai1[i1+k]), alpha, sig2, W)), 0, 1, 1, 1);
                        }
                        for (k = 0; k < ns0[i]; k++) {
                                alpha.at(d) = alpha_p;
                                logr += efi0[i0+k]*R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[i0+k]), alpha, sig2, W)), 0, 1, 0, 1);
                                alpha.at(d) = alpha_c;
                                logr -= efi0[i0+k]*R::pnorm(mu - abs(zeta(p, m, X.row(i), alpha, sig2, W) - zeta(p, m, X.row(sai0[i0+k]), alpha, sig2, W)), 0, 1, 0, 1);
                        }
                        i1 += ns1[i];
                        i0 += ns0[i];
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
}

// [[Rcpp::export]]
void  MCMC(const uvec& Y, const mat& X, const mat& W, const int& n, const int& p, const int& m, const int& n_sams, const int& n_burn, const int& n_skip, string prefix, string path_outs)
{
        // parameters
        double mu = 0.0;
        rowvec sig2 (p, fill::ones);
        rowvec alpha(m, fill::zeros);
        // MH parameters
        double del_mu = 1.0, mr_mu = 0.0;
        rowvec del_sig2 (p, fill::ones), mr_sig2 (p, fill::zeros);
        rowvec del_alpha(m, fill::ones), mr_alpha(m, fill::zeros);
        // open output files
        char* full;
        string nam;
        nam = "ll_chain"   ; full = mypaste0(path_outs, prefix, nam); ofstream ll_chain   ; ll_chain.open(full)   ;
        nam = "mu_chain"   ; full = mypaste0(path_outs, prefix, nam); ofstream mu_chain   ; mu_chain.open(full)   ;
        nam = "sig2_chain" ; full = mypaste0(path_outs, prefix, nam); ofstream sig2_chain ; sig2_chain.open(full) ;
        nam = "alpha_chain"; full = mypaste0(path_outs, prefix, nam); ofstream alpha_chain; alpha_chain.open(full);
        // chain
        int l, b, B = n_burn + n_skip*n_sams, n_disp = floor(0.1*double(B));
        for (b = 1; b <= B; b++) {
                // update
                sample_mu    (b, del_mu,    mr_mu,    n, p, m, mu, alpha, sig2, W, X, Y);
                sample_sig2  (b, del_sig2,  mr_sig2,  n, p, m, mu, alpha, sig2, W, X, Y);
                sample_alpha (b, del_alpha, mr_alpha, n, p, m, mu, alpha, sig2, W, X, Y);
                // save
                if ((b > n_burn) && (b % n_skip == 0)) {
                        ll_chain << loglik(n, p, m, mu, alpha, sig2, W, X, Y) << "\n";
                        mu_chain << mu << "\n";
                        for (l = 0; l < p; l++) sig2_chain  << sig2.at(l)  << "\n";
                        for (l = 0; l < m; l++) alpha_chain << alpha.at(l) << "\n";
                }
                // display
                if (b % n_disp == 0) {
                        Rcpp::Rcout << 100.0*double(b)/double(B) << "% completed" << endl;
                        // Rcpp::Rcout << 100.0*double(b)/double(B) << "% completed" << endl <<
                        // "Mixing rates" << endl <<
                        // " * mu    -> " << round_to_digits(mr_mu, 3)          << endl <<
                        // " * sig2  -> " << round_to_digits(mean(mr_sig2), 3)  << endl <<
                        // " * alpha -> " << round_to_digits(mean(mr_alpha), 3) << endl;
                }
        }
}

// [[Rcpp::export]]
void  MCMC_hat (const uvec& sai1, const uvec& sai0, const uvec& ns1, const uvec& ns0, const vec& efi0, const uvec& Y, const mat& X, const mat& W, const int& n, const int& p, const int& m, const int& n_sams, const int& n_burn, const int& n_skip, string prefix, string path_outs)
{
        // parameters
        double mu = 0.0;
        rowvec sig2 (p, fill::ones);
        rowvec alpha(m, fill::zeros);
        // MH parameters
        double del_mu = 1.0, mr_mu = 0.0;
        rowvec del_sig2 (p, fill::ones), mr_sig2 (p, fill::zeros);
        rowvec del_alpha(m, fill::ones), mr_alpha(m, fill::zeros);
        // open output files
        char* full;
        string nam;
        nam = "ll_chain"   ; full = mypaste0(path_outs, prefix, nam); ofstream ll_chain   ; ll_chain.open(full)   ;
        nam = "mu_chain"   ; full = mypaste0(path_outs, prefix, nam); ofstream mu_chain   ; mu_chain.open(full)   ;
        nam = "sig2_chain" ; full = mypaste0(path_outs, prefix, nam); ofstream sig2_chain ; sig2_chain.open(full) ;
        nam = "alpha_chain"; full = mypaste0(path_outs, prefix, nam); ofstream alpha_chain; alpha_chain.open(full);
        // chain
        int l, b, B = n_burn + n_skip*n_sams, n_disp = floor(0.1*double(B));
        for (b = 1; b <= B; b++) {
                // update
                sample_mu_hat    (sai1, sai0, ns1, ns0, efi0, b, del_mu,    mr_mu,    n, p, m, mu, alpha, sig2, W, X, Y);
                sample_sig2_hat  (sai1, sai0, ns1, ns0, efi0, b, del_sig2,  mr_sig2,  n, p, m, mu, alpha, sig2, W, X, Y);
                sample_alpha_hat (sai1, sai0, ns1, ns0, efi0, b, del_alpha, mr_alpha, n, p, m, mu, alpha, sig2, W, X, Y);
                // save
                if ((b > n_burn) && (b % n_skip == 0)) {
                        ll_chain << loglik_hat(sai1, sai0, ns1, ns0, efi0, n, p, m, mu, alpha, sig2, W, X) << "\n";
                        mu_chain << mu << "\n";
                        for (l = 0; l < p; l++) sig2_chain  << sig2.at(l)  << "\n";
                        for (l = 0; l < m; l++) alpha_chain << alpha.at(l) << "\n";
                }
                // display
                if (b % n_disp == 0) {
                        Rcpp::Rcout << 100.0*double(b)/double(B) << "% completed" << endl;
                        // Rcpp::Rcout << 100.0*double(b)/double(B) << "% completed" << endl <<
                        // "Mixing rates" << endl <<
                        // " * mu    -> " << round_to_digits(mr_mu, 3)          << endl <<
                        // " * sig2  -> " << round_to_digits(mean(mr_sig2), 3)  << endl <<
                        // " * alpha -> " << round_to_digits(mean(mr_alpha), 3) << endl;
                }
        }
}