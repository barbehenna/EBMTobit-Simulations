# Simulations Comparing Left-Censored Imputation Methods

## What is this?

This repository contains simulation scripts and artifacts accompanying Barbehenn and Zhao (2023) "Nonparametric empirical Bayes biomarker imputation and estimation" (https://arxiv.org/abs/2306.07239). The following methods are compared:

- EBM-Tobit [^1]
- Tobit maximum likelihood estimator [^1]
- GSimp [^2]
- QRILC [^3]
- zCompositions [^4]
- trKNN [^5]
- Empirical Bayes g-model with generalized exemplar support [^1]
- Empirical Bayes g-model with oracle support [^1]
- Vectorized empirical Bayes g-model with oracle support [^1]

The last two methods cannot be used in practice but are used to motivate the empirical Bayes framework for this imputation problem. 


## How to Set-up

Follow the following steps to set up the simulation environment. We used R version 4.2.2 (2022-10-31) in our simulations.

1. Set up the GSimp method
    a. Clone the https://github.com/WandeRum/GSimp repository
    a. Update line 3 in this repository's `Utils.R` to contain the location of the GSimp repository
    a. Remove the following lines: lines 9-10 in `GSimp.R`, lines 7-8 in `GSimp_evaluation.R`, and lines 7-8 in `Impute_wrapper.R`
1. Ensure that the following packages are installed (later versions may work, but the following versions were used in the published simulations); most of the packages are on CRAN, but some are on BioConductor
```
 [1] ropls_1.30.0          reshape2_1.4.4        vegan_2.6-4          
 [4] lattice_0.21-8        permute_0.9-7         doParallel_1.0.17    
 [7] iterators_1.0.14      foreach_1.5.2         imputeLCMD_2.1       
[10] impute_1.72.3         pcaMethods_1.90.0     Biobase_2.58.0       
[13] BiocGenerics_0.44.0   norm_1.0-11.0         tmvtnorm_1.5         
[16] gmm_1.7               sandwich_3.0-2        mvtnorm_1.1-3        
[19] abind_1.4-5           missForest_1.5        magrittr_2.0.3       
[22] FNN_1.1.3.2           rpart_4.1.19          glmnet_4.1-7         
[25] Matrix_1.5-4.1        randomForest_4.7-1.1  ggplot2_3.4.2        
[28] data.table_1.14.8     zCompositions_1.4.0-1 truncnorm_1.0-9      
[31] NADA_1.6-1.1          survival_3.5-5        MASS_7.3-60          
[34] ebTobit_1.0.1        
```


## How to Run Simulations

After setting up all of the simulations, use
```
Rscript Run_Simulations.R
```
to run the simulations. The number of simulations and combinations of censoring levels are defined on lines 122-128 of `Run_Simulations.R`. An additional log file containing environmental information is produced before simulations start.


## Files

- `EMBTobit.R`: an R function implementing Algorithm 1 (the EBM-Tobit estimator in https://arxiv.org/abs/2306.07239)
- `Utils.R`: a helper R script to load the GSimp method and generate data for simulations
- `Run_Simulations.R`: an R script to reproducibly set up and run the simulations; both `sims.csv` and `sessionInfo.log` are produced by this script
- `sims.csv`: simulation results used for method comparison 
- `sessionInfo.log`: a text file containing the full output of `sessionInfo()` used in for the simulations


## References

[^1]: Barbehenn, A., & Zhao, S. D. (2023). Nonparametric empirical Bayes biomarker imputation and estimation.
[^2]: Wei R, Wang J, Jia E, Chen T, Ni Y, et al. (2018) GSimp: A Gibbs sampler based left-censored missing value imputation approach for metabolomics studies. PLOS Computational Biology 14(1): e1005973. https://doi.org/10.1371/journal.pcbi.1005973
[^3]: Lazar, C., & Burger, T. (2022). imputeLCMD: A Collection of Methods for Left-Censored Missing Data Imputation. https://CRAN.R-project.org/package=imputeLCMD
[^4]: Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015). zCompositions — R package for multivariate imputation of left-censored data under a compositional approach. Chemometrics and Intelligent Laboratory Systems, 143, 85–96. https://doi.org/https://doi.org/10.1016/j.chemolab.2015.02.019
[^5]: Shah, J. S., Rai, S. N., DeFilippis, A. P., Hill, B. G., Bhatnagar, A., & Brock, G. N. (2017). Distribution based nearest neighbor imputation for truncated high dimensional data with applications to pre-clinical and clinical metabolomics studies. BMC Bioinformatics, 18(1), 114. https://doi.org/10.1186/s12859-017-1547-6
