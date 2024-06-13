#### Load methods and data ####

library(ebTobit)
library(zCompositions)
library(data.table)
source("Utils.R") # GSimp and `generate_data()`
source("EBMTobit.R")
source("fullEM.R")

TASK_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("TASK_ID =", TASK_ID))

#### Simulation Helpers ####

#' Simulation Runner: fit data with all methods
#' @param L n x p matrix of lower bounds; L <= X <= R
#' @param R n x p matrix of upper bounds
#' @param TH n x p matrix of true means; X ~ p(x | t = TH)
#' @param Xobs n x p matrix of uncensored observations
fit_models = function(L, R, TH, Xobs) {
    # helpers
    n = nrow(L)
    p = ncol(L)
    X = L; X[L < R] = NA

    # Collect method fits
    est = list()
	
	
	#### empirical Bayes oracles ####

    ## EB oracle posterior mean
    est[["eb_oracle_knots"]] = fitted(ebTobit(L = L, R = R, gr = TH, algorithm = ALG, rtol = RTOL))

    ## EB Vectorized oracle
    if (n <= 1000) {
        est[["eb_vectorized_oracle"]] = matrix(fitted(ebTobit(L = c(L), R = c(R), gr = c(TH))), n, p)
    }

	
	#### empirical Bayes methods ####
	
    # MCMC sample from marginal
    est[["eb_EBMTobit"]] = EBMTobit(L = L, R = R, K = 50, algorithm = ALG, rtol = RTOL)
	
    # generalized exemplar
    est[["eb_exemplar_MLE"]] = fitted(ebTobit(L = L, R = R, gr = (L + R)/2, algorithm = ALG, rtol = RTOL))
    
    # Full EM
    est[["eb_full_mle_seed"]] = fitted(fullEM(L = L, R = R, gr = (L + R)/2, algorithm = ALG, rtol = RTOL))

    # full EM fitting method
    mods = lapply(1:50, FUN = function(i) fullEM(L = L, R = R, gr = sapply(1:ncol(L), function(j) runif(200, min = min(L[,j]), max = max(R[,j]))), rtol = RTOL))
    est[["eb_fullSingle"]] = fitted(mods[[1]])
    est[["eb_fullBest"]] = fitted(mods[[which.max(sapply(mods, logLik))]])
    est[["eb_fullAvg50"]] = apply(simplify2array(lapply(mods, fitted)), MARGIN = 1:2, FUN = mean)
	
    #### other methods ####
    
    if (!all(L == R)) {
        ## Tobit MLE (full rank)
        est[["TobitMLE"]] = (L + R)/2
        
        ## QRILC
        est[["QRILC"]] = impute.QRILC(X)[[1]]
        
        ## zCompositions
        LOD = exp(apply(R, 2, min, na.rm = TRUE))
        LOD[colSums(is.na(X)) == 0] = 0
        est[["zCompositions"]] = log(as.matrix(multLN(X = exp(X), label = NA, dl = LOD)))
        
        ## HM
        est[["HalfMin"]] = apply(X, 2, function(x) {
            x[is.na(x)] <- min(x, na.rm = T)/2
            x
        })
        
        ## trKNN
        est[["trKNN"]] = tryCatch(
            t(imputeKNN(t(X), k = 3, dist = "truncation", perc = 0)),
            error = function(e) matrix(NA, n, p)
        )
        
        ## SVD
        full_mod = pca(X, method = "svdImpute", nPcs = p, cv = "q2")@cvstat
        est[["SVD"]] = pca(X, method = "svdImpute", nPcs = which.min(full_mod))@completeObs
        
        ## GSimp
        est[["GSimp"]] = GS_impute(as.data.frame(X), iters_each=50, iters_all=10,
                                   initial='qrilc', lo=apply(L, 2, min) - 0.1, hi='min', n_cores=1,
                                   imp_model='glmnet_pred')$data_imp |> as.matrix()
    }
    
    # collect results into a data.table
    df = do.call("rbind", lapply(names(est), function(method.name) {
        fit = est[[method.name]]
        data.table(
            method = method.name,
            expand.grid(
                "using" = c("all", "cen"),
                "metric" = c("mse", "nrmse", "pearson", "spearman"),
                "target" = c("TH", "X")
            ),
            value  = c(
                mean((TH - fit)^2),
                mean((TH[is.na(X)] - fit[is.na(X)])^2),
                sqrt(mean((TH - fit)^2)/var(c(TH))),
                missForest::nrmse(ximp = fit, xmis = X, xtrue = TH),
                cor(c(TH), c(fit), method = "pearson"),
                cor(c(TH[is.na(X)]), c(fit[is.na(X)]), method = "pearson"),
                cor(c(TH), c(fit), method = "spearman"),
                cor(c(TH[is.na(X)]), c(fit[is.na(X)]), method = "spearman"),
                mean((Xobs - fit)^2),
                mean((Xobs[is.na(X)] - fit[is.na(X)])^2),
                sqrt(mean((Xobs - fit)^2)/var(c(TH))),
                missForest::nrmse(ximp = fit, xmis = X, xtrue = Xobs),
                cor(c(Xobs), c(fit), method = "pearson"),
                cor(c(Xobs[is.na(X)]), c(fit[is.na(X)]), method = "pearson"),
                cor(c(Xobs), c(fit), method = "spearman"),
                cor(c(Xobs[is.na(X)]), c(fit[is.na(X)]), method = "spearman")
            )
        )
    }))
	df
}

run.sim = function(n, p, mc, mp) {
    d = generate_data(n = n, p = p, mc = mc, mp = mp)
    fit_models(L = d$L, R = d$R, TH = d$TH, Xobs = d$Xobs)
}


#### Run Simulations ####

# Reproducibility: record simulation session information
writeLines(capture.output(sessionInfo()), "sessionInfo.log")

# Reproducibility: set seeds and threading
set.seed(TASK_ID)
setDTthreads(TASK_ID)

sims = as.data.table(expand.grid(
    iter = TASK_ID,
    n = c(100, 1000),
    p = c(10, 25),
    mc = c(0.1, 0.3, 0.5),
    mp = c(0.1, 0.3, 0.5)
))

# Run simulations
sims = sims[, run.sim(n = n, p = p, mc = mc, mp = mp), by = .(iter, n, p, mc, mp)]
fwrite(sims, paste0("sims_", TASK_ID, ".csv"))
summary(warnings())


# Run supplemental simulations - noise distribution
run.sim.noise = function(n, p, noise, mc, mp) {
    d = generate_data_noise(n = n, p = p, noise = noise, mc = mc, mp = mp)
    fit_models(L = d$L, R = d$R, TH = d$TH, Xobs = d$Xobs)
}

sims.noise = as.data.table(expand.grid(
    iter = TASK_ID,
    n = c(1000),
    p = c(10,25),
    noise = c("gaussian", "exponential", "uniform"),
    mc = c(0.3),
    mp = c(0.1)
))

sims.noise = sims.noise[, run.sim.noise(n = n, p = p, noise = as.character(noise), mc = mc, mp = mp), by = .(iter, n, p, noise, mc, mp)]
fwrite(sims.noise, paste0("sims_noise_", TASK_ID, ".csv"))
summary(warnings())


# Run supplemental simulations - rank (and scale to some extent since theta is chi-squared k)
run.sim.rank = function(n, p, k, mc, mp) {
    d = generate_data_rank(n = n, p = p, k = k, mc = mc, mp = mp)
    fit_models(L = d$L, R = d$R, TH = d$TH, Xobs = d$Xobs)
}

sims.rank = as.data.table(expand.grid(
    iter = TASK_ID,
    n = c(1000),
    p = c(25),
    k = c(1, 5, 10, 20),
    mc = c(0.3),
    mp = c(0.1)
))

sims.rank = sims.rank[, run.sim.rank(n = n, p = p, k = k, mc = mc, mp = mp), by = .(iter, n, p, k, mc, mp)]
fwrite(sims.rank, paste0("sims_rank_", TASK_ID, ".csv"))
summary(warnings())


# Run supplemental simulations - uncensored algorithm comparison
run.sim.alg = function(n, p) {
    TH = matrix(rnorm(n*p), n, p)
    Xobs = TH + matrix(rnorm(n*p), n, p)
    
    fit_models(L = Xobs, R = Xobs, TH = TH, Xobs = Xobs)
}

sims.alg = as.data.table(expand.grid(
    iter = TASK_ID,
    n = c(100,1000),
    p = c(5,10,25)
))

sims.alg = sims.alg[, run.sim.alg(n = n, p = p), by = .(iter, n, p)]
fwrite(sims.alg, paste0("sims_alg_", TASK_ID, ".csv"))
summary(warnings())

print(paste("Finished TASK_ID", TASK_ID))
