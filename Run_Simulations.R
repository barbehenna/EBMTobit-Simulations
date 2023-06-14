#### Load methods and data ####

library(ebTobit)
library(zCompositions)
library(data.table)
library(ggplot2)
source("Utils.R") # GSimp and `generate_data()`
source("EBMTobit.R")


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
    est[["eb_vectorized_oracle"]] = matrix(fitted(ebTobit(L = c(L), R = c(R), gr = c(TH))), n, p)

	
	#### empirical Bayes methods ####
	
    # MCMC sample from marginal
    est[["eb_EBMTobit"]] = EBMTobit(L = L, R = R, K = 50, algorithm = ALG, rtol = RTOL)
	
    # generalized exemplar
    est[["eb_exemplar_MLE"]] = fitted(ebTobit(L = L, R = R, gr = (L + R)/2, algorithm = ALG, rtol = RTOL))

	
	#### other methods ####
	
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

    ## GSimp
    est[["GSimp"]] = GS_impute(as.data.frame(X), iters_each=50, iters_all=10,
                               initial='qrilc', lo=apply(L, 2, min) - 0.1, hi='min', n_cores=1,
                               imp_model='glmnet_pred')$data_imp |> as.matrix()



    # collect results into a data.table
    df = do.call("rbind", lapply(names(est), function(method.name) {
        fit = est[[method.name]]
        data.table(
            method = method.name,
            expand.grid(
                "using" = c("all", "cen"),
                "metric" = c("mse", "pearson", "spearman"),
                "target" = c("TH", "X")
            ),
            value  = c(
                mean((TH - fit)^2),
                mean((TH[is.na(X)] - fit[is.na(X)])^2),
                cor(c(TH), c(fit), method = "pearson"),
                cor(c(TH[is.na(X)]), c(fit[is.na(X)]), method = "pearson"),
                cor(c(TH), c(fit), method = "spearman"),
                cor(c(TH[is.na(X)]), c(fit[is.na(X)]), method = "spearman"),
                mean((Xobs - fit)^2),
                mean((Xobs[is.na(X)] - fit[is.na(X)])^2),
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
set.seed(1)
setDTthreads(1)

sims = as.data.table(expand.grid(
    iter = 1:200,
    n = 1000,
    p = 25,
    mc = c(0.1, 0.3, 0.5),
    mp = c(0.1, 0.3, 0.5)
))

# Run simulations
sims = sims[, run.sim(n = n, p = p, mc = mc, mp = mp), by = .(iter, n, p, mc, mp)]
fwrite(sims, "sims.csv")
summary(warnings())
