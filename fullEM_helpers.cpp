// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include <cmath>

// ---- Helper Functions ----
// Compute the centered Gaussian pdf
inline double phi(double x, double s = 1.0) {
    return exp(-0.5 * x * x / s / s) * M_SQRT1_2 * M_2_SQRTPI / 2 / s;
}

// Compute the *difference* in standard Gaussian cdfs (m = 0, s = 1)
// Phi(x) - Phi(y)
inline double dPhi(double x, double y) {
    return (erf(x * M_SQRT1_2) - erf(y * M_SQRT1_2)) / 2;
}


// ---- Conditional Mean Calculation ----

// Compute the mean of a truncated normal random variable
inline double tr_mean(double m, double s, double a, double b) {
    double alpha = (a - m) / s;
    double beta  = (b - m) / s;
    
    // special case: E[X | a <= X <= b=a] = E[X | X=a] = a
    if (a == b) {
        return a;
    }
    
    // numerical stability heuristics (dPhi is approximately 0)
    if (alpha > 6) {
        return a;
    }
    
    if (beta < -6) {
        return b;
    }
    
    // mean calculation
    return m - s * (phi(beta) - phi(alpha)) / dPhi(beta, alpha);
}

// RcppParallel worker -- compute the conditional mean
// Parallelize the calculation by row of L (sample/observation)
// This worker fills a vector of length n*p*m so that it can be later converted to an array.
struct CondMean : public RcppParallel::Worker {
     
    // Input matrices
    const RcppParallel::RMatrix<double> L;
    const RcppParallel::RMatrix<double> R;
    const RcppParallel::RMatrix<double> G;
    const RcppParallel::RMatrix<double> s1;

    // Output matrix
    RcppParallel::RVector<double> cmean;

    // Initialize from Rcpp input and output matrices
    CondMean(const Rcpp::NumericMatrix L, const Rcpp::NumericMatrix R,
             const Rcpp::NumericMatrix G, const Rcpp::NumericMatrix s1,
             const Rcpp::NumericVector cmean)
        : L(L), R(R), G(G), s1(s1), cmean(cmean) {}

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t j = 0; j < L.ncol(); j++) {
                for (std::size_t k = 0; k < G.nrow(); k++) {
                    cmean[i + j*L.nrow() + k*L.nrow()*L.ncol()] = tr_mean(G(k,j), s1(i,j), L(i,j), R(i,j));
                }
            }
        }
    }
     
};

//' Helper Function - generate conditional mean array
//'
//' Compute an array E whose entries are 
//' E[i,j,k] = E[X_ij | gr = k, L_ij <= X_ij <= R_ij ] for
//' observed interval (L_ij, R_ij) on measurement X_ij and grid of means gr_kj.
//'
//' @param L n x p matrix of lower bounds
//' @param R n x p matrix of upper bounds
//' @param gr m x p matrix of candidate means
//' @param s1 n x p matrix of standard deviations
//' @return the n x p x m expectation array under partial interval censoring
//'
//' @examples
//' # set-up
//' n = 100; m = 50; p = 5
//' gr = matrix(stats::rnorm(m*p), m, p)
//' L = R = matrix(stats::rnorm(n*p), n, p)
//' s1 = matrix(1, n, p)
//' missing.idx = sample.int(n = n*p, size = p*p)
//' L[missing.idx] = L[missing.idx] - stats::runif(p, 0, 1)
//'
//' # R solution
//' y = array(dim = c(n, p, m))
//' for (i in 1:n) {
//'     for (j in 1:p) {
//'         for (k in 1:m) {
//'             if (L[i,j] == R[i,j]) y[i,j,k] = L[i,j]
//'             else if (L[i,j] - gr[k,j] > 6) y[i,j,k] = L[i,j]
//'             else if (R[i,j] - gr[k,j] < -6) y[i,j,k] = R[i,j]
//'             else y[i,j,k] = gr[k,j] - (dnorm(R[i,j] - gr[k,j]) - dnorm(L[i,j] - gr[k,j])) / (pnorm(R[i,j] - gr[k,j]) - pnorm(L[i,j] - gr[k,j]))
//'         }
//'     }
//' }
//'
//' # Compare R to RcppParallel method
//' all.equal(y, ebTobit:::compMeans(L, R, gr, s1))
//' @useDynLib ebTobit
//' @importFrom Rcpp evalCpp
//' @import RcppParallel
// [[Rcpp::export]]
Rcpp::NumericVector compMeans(Rcpp::NumericMatrix L, 
                              Rcpp::NumericMatrix R, 
                              Rcpp::NumericMatrix gr, 
                              Rcpp::NumericMatrix s1) {

    // allocate the vector we will return and record dimensions of output array
    Rcpp::NumericVector out(L.nrow()*L.ncol()*gr.nrow());
    Rcpp::NumericVector dim(3);   
    dim[0] = L.nrow();
    dim[1] = L.ncol();
    dim[2] = gr.nrow();
    
    // create the worker
    CondMean condmean(L, R, gr, s1, out);

    // call worker with parallelFor
    // parallel over rows in out
    RcppParallel::parallelFor(0, L.nrow(), condmean);
    
    // convert out to array
    out.attr("dim") = dim;
    
    return out;
}


// ---- Conditional Variance Calculation ----

// Compute the variance of a truncated normal random variable
inline double tr_var(double m, double s, double a, double b) {
    double alpha = (a - m) / s;
    double beta  = (b - m) / s;
    
    // special case: E[X^2 | a <= X <= b=a] = E[X^2 | X=a] = a^2
    if (a == b) {
        return a*a;
    }
    
    // numerical stability heuristics (dPhi is approximately 0)
    if (alpha > 6) {
        return a*a + s*s;
    }
    
    if (beta < -6) {
        return b*b + s*s;
    }
    
    // variance calculation
    return s*s*(1.0 - (beta * phi(beta) - alpha * phi(alpha))/dPhi(beta, alpha) - pow((phi(beta) - phi(alpha))/dPhi(beta, alpha), 2));
}

// RcppParallel worker -- compute the conditional variance
// Parallelize the calculation by row of L (sample/observation)
// This worker fills an n x m matrix
struct CondVar : public RcppParallel::Worker {
    
    // Input matrices
    const RcppParallel::RMatrix<double> L;
    const RcppParallel::RMatrix<double> R;
    const RcppParallel::RMatrix<double> G;
    const RcppParallel::RMatrix<double> s1;
    
    // Output matrix
    RcppParallel::RMatrix<double> cvar;
    
    // Initialize from Rcpp input and output matrices
    CondVar(const Rcpp::NumericMatrix L, const Rcpp::NumericMatrix R,
            const Rcpp::NumericMatrix G, const Rcpp::NumericMatrix s1,
            const Rcpp::NumericMatrix cvar)
        : L(L), R(R), G(G), s1(s1), cvar(cvar) {}
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t k = 0; k < G.nrow(); k++) {
                double tmp = 0;
                for (std::size_t j = 0; j < G.ncol(); j++) {
                    tmp += tr_var(G(k,j), s1(i,j), L(i,j), R(i,j));
                }
                cvar(i,k) = tmp;
            }
        }
    }
    
};


//' Helper Function - generate conditional variance matrix
//'
//' Compute a matrix V whose entries are
//' V[i,k] = var[X_i* | gr = k, L_i* <= X_i* <= R_i* ] for
//' observed interval (L_i*, R_i*) on measurement X_i* and grid of means gr_k*.
//'
//' @param L n x p matrix of lower bounds
//' @param R n x p matrix of upper bounds
//' @param gr m x p matrix of candidate means
//' @param s1 n x p matrix of standard deviations
//' @return the n x p x m expectation array under partial interval censoring
//'
//' @examples
//' # set-up
//' n = 100; m = 50; p = 5
//' gr = matrix(stats::rnorm(m*p), m, p)
//' L = R = matrix(stats::rnorm(n*p), n, p)
//' s1 = matrix(1, n, p)
//' missing.idx = sample.int(n = n*p, size = p*p)
//' L[missing.idx] = L[missing.idx] - stats::runif(p, 0, 1)
//'
//' # R solution
//' v = array(dim = c(n, p, m))
//' for (i in 1:n) {
//'     for (j in 1:p) {
//'         for (k in 1:m) {
//'             if (L[i,j] == R[i,j]) v[i,j,k] = L[i,j]^2
//'             else if (L[i,j] - gr[k,j] > 6) v[i,j,k] = L[i,j]^2 + 1
//'             else if (R[i,j] - gr[k,j] < -6) v[i,j,k] = R[i,j]^2 + 1
//'             else {
//'                 l = L[i,j] - gr[k,j]
//'                 r = R[i,j] - gr[k,j]
//'                 v[i,j,k] =  1 - (r*dnorm(r) - l*dnorm(l))/(pnorm(r) - pnorm(l)) - ((dnorm(r) - dnorm(l))/(pnorm(r) - pnorm(l)))^2
//'             }
//'         }
//'     }
//' }
//' # Compare R to RcppParallel method
//' all.equal(apply(v, c(1,3), sum), ebTobit:::compVars(L, R, gr, s1))
//' @useDynLib ebTobit
//' @importFrom Rcpp evalCpp
//' @import RcppParallel
// [[Rcpp::export]]
Rcpp::NumericVector compVars(Rcpp::NumericMatrix L,
                             Rcpp::NumericMatrix R,
                             Rcpp::NumericMatrix gr,
                             Rcpp::NumericMatrix s1) {

    // allocate the matrix we will return
    Rcpp::NumericMatrix out(L.nrow(), gr.nrow());

    // create the worker
    CondVar condvar(L, R, gr, s1, out);

    // call worker with parallelFor
    // parallel over rows in out
    RcppParallel::parallelFor(0, L.nrow(), condvar);

    return out;
}
