#### Pool csv files from Slurm Array Job ####

# Each set of csv files are merged into one large one
# the originals are copied into a tar ball then deleted

library(data.table)


## Main simulations
fnames = list.files(pattern = "sims_\\d+.csv")
fwrite(rbindlist(lapply(fnames, fread)), file = "sims.csv")
tar("sims.tar.gz", files = fnames, compression = "gzip")
file.remove(fnames)


## Noise distribution simulations
fnames = list.files(pattern = "sims_noise_\\d+.csv")
fwrite(rbindlist(lapply(fnames, fread)), file = "sims_noise.csv")
tar("sims_noise.tar.gz", files = fnames, compression = "gzip")
file.remove(fnames)


## structure simulations
fnames = list.files(pattern = "sims_rank_\\d+.csv")
fwrite(rbindlist(lapply(fnames, fread)), file = "sims_rank.csv")
tar("sims_rank.tar.gz", files = fnames, compression = "gzip")
file.remove(fnames)


## algorithm simulations
fnames = list.files(pattern = "sims_alg_\\d+.csv")
fwrite(rbindlist(lapply(fnames, fread)), file = "sims_alg.csv")
tar("sims_alg.tar.gz", files = fnames, compression = "gzip")
file.remove(fnames)
