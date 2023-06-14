#### Set-Up and Helper Functions ####

GSIMP = "~/git/GSimp" # location of https://github.com/WandeRum/GSimp
                      # (commit: e2a5d56fb33846f139635fe0209a0c374694564e)
                      # modified to not source other scripts in the directory. Remove the following lines:
					  #   - `GSimp.R` lines 9-10
					  #   - `GSimp_evaluation.R` lines 7-8
					  #   - `Impute_wrapper.R` lines 7-8
ALG = "EM"
RTOL = 1e-6

base_data = log(as.matrix(ebTobit::BileAcid))

source(file.path(GSIMP, "Trunc_KNN/Imput_funcs.r"))
source(file.path(GSIMP, "Prediction_funcs.R"))
source(file.path(GSIMP, "MVI_global.R"))
source(file.path(GSIMP, "GSimp.R"))
source(file.path(GSIMP, "Impute_wrapper.R"))
source(file.path(GSIMP, "GSimp_evaluation.R"))


#' Generate data with independent censoring
#' Since X is random, it's possible there's censoring even when mp = mc = 0.
#' @param setting character denoting structure of data; either "realcov" or "rank*"
#' (where * is taken to be the latent rank)
#' @param n number of samples (number of rows)
#' @param p dimension of samples (number of columns)
#' @param mc percent of columns that have missing values (randomly sample round(mc*p) columns)
#' @param mp lower oracle quantile used to threshold
generate_data = function(n = 1000, p = 10, mc = 0.1, mp = 0.1) {
    print(paste("---- setting --", n, p, mc, mp, "----"))
    
    # generate true means
    cols = sample(1:ncol(base_data), size = p, replace = TRUE)
    MU = colMeans(base_data[,cols])
    eig_ = eigen(cov(base_data[,cols]))
    L = eig_$vectors %*% diag(sqrt(pmax(eig_$values, 0))) # all.equal(unname(cov(base_data[,cols])), tcrossprod(L))
    TH = t(L %*% matrix(rnorm(n*p), p, n) + MU)
    
    # generate noisy observations
    X = TH + matrix(rnorm(n*p), n, p)
    
    # add left censoring
    cen = sample(1:p, size = round(mc*p)) # columns with censoring
    print(paste("censoring columns:", paste(cen, collapse = ", ")))
    cen_bdd = apply(TH, 2, quantile, prob = mp)
    cen_bdd[-cen] = -Inf
    L = sapply(1:p, function(j) ifelse(X[,j] < cen_bdd[j], min(TH[,j]) - 6*sd(TH[,j]), X[,j]))
    R = sapply(1:p, function(j) ifelse(X[,j] < cen_bdd[j], cen_bdd[j], X[,j]))
    
    list(L = L, R = R, TH = TH, Xobs = X)
}
