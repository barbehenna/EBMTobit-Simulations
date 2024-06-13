library(ebTobit)
Rcpp::sourceCpp("fullEM_helpers.cpp")

#' @rdname ebTobit
#' @param maxiter early stopping condition
#' @param rtol convergence criteria: abs(loss_new - loss_old)/abs(loss_new)
#' @export
fullEM <- function(L, R = L, gr = (R+L)/2, s1 = 1, pos_lik = TRUE, 
                         maxiter = 10000L, rtol = 1e-6, ...) {
    
    # expand s1 to match L
    if (length(s1) == 1) {
        s1 <- matrix(s1, nrow = nrow(L), ncol = ncol(L))
    }
    
    # basic checks
    stopifnot(is.matrix(L))
    stopifnot(all(dim(L) == dim(R)))
    stopifnot(all(L <= R))
    stopifnot(is.matrix(gr))
    stopifnot(ncol(gr) == ncol(L))
    stopifnot(all(dim(s1) == dim(L)))
    stopifnot(all(s1 > 0))
    
    # set-up
    pr <- rep(1, nrow(gr))/nrow(gr)
    loss_old <- -.Machine$double.xmax
    
    # Full EM fitting (gr and pr)
    for (iter in seq_len(maxiter)) {
        ## E step
        lik <- likMat(L, R, gr, s1)
        if (pos_lik) lik <- pmax(lik, .Machine$double.xmin, na.rm = TRUE)
        
        z <- lik %*% diag(pr) # conditional probabilities
        z <- z / rowSums(z)
        y = compMeans(L, R, gr, s1) # conditional means
        v = compVars(L, R, gr, s1) # conditional variances

        ## M step
        pr <- colMeans(z)
        for (k in seq_along(pr)) {
            gr[k,] <- colSums(y[,,k] * z[,k], na.rm = TRUE) / sum(z[,k], na.rm = TRUE)
        }
        loss_new <- nrow(L)*sum(pr*log(pr)) - sum(v %*% pr)/2
        
        ## Convergence check (weights stop updating)
        if ( is.na(loss_new) | abs((loss_new - loss_old)/loss_new) < rtol ) break
        loss_old <- loss_new
    }
    
    lik = likMat(L, R, gr, s1)
    if (pos_lik) lik = pmax(lik, .Machine$double.xmin, na.rm = TRUE)

    new_ebTobit(
        prior = pr,
        gr = gr,
        lik = lik
    )
}
