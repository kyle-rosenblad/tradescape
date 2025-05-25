buildtradescape <- function(allele_count, # vector of integer values; each element is the reference allele count for one individual (from one site) at one snp
                            snp_idx, # integer vector of same length as allele_count; index variable to indicate which individual each allele count pertains to
                            site_idx, # integer vector of same length as allele_count; index variable to indicate which site each allele count pertains to
                            envdata, # matrix of real values; rows are individuals (NOT SITES), columns are environmental variables. maybe make rows sites in future release since this might be more intuitive for users.
                            coords, # matrix of real values; rows are sites, columns are x and y spatial coordinates
                            ploidy, # single integer value; defines max allele count per individual
                            cores=2, # how many processor cores to use in parallel for sampling
                            chains=2, # how many Hamiltonian Monte Carlo (HMC) sampler chains to run
                            iter=2000, # how many total HMC samples to draw
                            warmup=1000, # how many HMC samples to use as warmup
                            adapt_delta=0.9 # parameter that adjusts how "cautiously" the HMC sampler explores the posterior geometry. Greater values (closer to 1) mean more caution.
                            ){
  # for each snp, flip a virtual coin to decide whether to switch the (arbitrary) reference allele designation
  unique_snps <- unique(snp_idx)
  flip <- sample(c(FALSE, TRUE), length(unique_snps), replace = TRUE)
  snps_to_flip <- flip[match(snp_idx, unique_snps)]
  allele_count[snps_to_flip] <- ploidy - allele_count[snps_to_flip]
  
  # define some additional needed inputs for the Stan code
  N_obs <- length(allele_count)
  N_snps <- length(unique(snp_idx))
  N_sites <- length(unique(site_idx))
  N_envars <- ncol(envdata)
  
  # put necessary inputs in list format for Stan
  stan_data <- list(N_snps,
                    N_obs,
                    N_sites,
                    N_envars,
                    snp_idx,
                    site_idx,
                    allele_count,
                    envdata,
                    coords)

  stanmod <- rstan::sampling(stanmodels$tradescape_binomial,
                             data=stan_data,
                             cores=N_cores,
                             chains=N_chains,
                             warmup=warmup,
                             iter=iter,
                             verbose=F,
                             control=list(adapt_delta=adapt_delta))
  return(stanmod)
}

