#' Empirical Bayes matrix estimation under a Tobit likelihood
#' 
#' A heuristic method for estimating the posterior mean in nonparametric
#' g-modeling problems. Initial support is taken from the user (by default the
#' maximum likelihood estimates are used), on consecutive iterations samples drawn
#' from an approximate marginal distribution are used as support for the g
#' prior. It is assumed throughout this method that the noise has unit standard
#' deviation.
#' 
#' @param L n x p matrix of lower bounds
#' @param R n x p matrix of upper bounds
#' @param gr.init m x p matrix containing m support points for the p-dimensional prior
#' @param K number of iterations
#' @param ... additional parameters passed to \code{ebTobit::ebTobit} such as rtol
#' @returns n x p matrix of estimated means
#' @import ebTobit
#' @examples 
#' set.seed(1)
#' n = 1000; p = 3
#' TH = matrix(stats::rexp(n*p, rate = 0.2), n, p)
#' X = TH + matrix(stats::rnorm(n*p), n, p)
#' ldl = 0.2 # observe [0, ldl] if X falls below ldl
#' L = ifelse(X < ldl, 0, X)
#' R = ifelse(X < ldl, ldl, X)
#' est = EBMTobit(L, R)
#' mean((est[L<R] - TH[L<R])^2) # imputation
EBMTobit <- function(L, R = L, gr.init = (L+R)/2, K = 50, ...) {
    # preliminary checks
    stopifnot(all(dim(L) == dim(R)))
    stopifnot(ncol(L) == ncol(gr.init))

    # set-up
    n <- nrow(L)
    p <- ncol(L)
    m <- nrow(gr.init)
    out <- matrix(0, n, p)
    gr.next <- gr.init

    # iterative fitting
    for (i in seq_len(K)) {
        fit <- ebTobit(L = L, R = R, gr = gr.next, ...)
        gr.next <- fit$gr[sample(m, size = m, replace = TRUE, prob = fit$prior),] + matrix(rnorm(m*p), m, p)
        out <- out + fitted(fit)
    }
    
    out / K
}
