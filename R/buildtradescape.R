#' Hierarchical Bayesian modeling of local adaptation to multidimensional environmental gradients
#'
#' @export
#' @param allele_count Integer vector, in which each element represents the count of the reference allele for one individual at one SNP.
#' @param snp_idx Integer vector of same length as allele_count, serving as an index that indicates which SNP each value of allele_count pertains to.
#' @param site_idx Integer vector of same length as allele_count, serving as an index that indicates which geographic sampling site each value of allele_count pertains to.
#' @param envdata Matrix of real values, in which rows represent individual organisms--not sampling sites--and columns represent environmental variables to be used as predictors of allele frequencies.
#' @param coords Matrix of real values, in which rows represent geographic sampling sites, and columns represent x and y geographic coordinates.
#' @param ploidy Single integer value representing the ploidy of the individuals in the data set--i.e., the maximum possible value for each element in allele_count.
#' @param cores Number of processor cores to use for parallel Hamiltonian Monte Carlo (HMC) posterio sampling.
#' @param chains Number of sampler chains for Stan's Hamiltonian Monte Carlo (HMC) sampler.
#' @param warmup Number of warmup iterations for Stan's Hamiltonian Monte Carlo (HMC) sampler.
#' @param iter Total number of posterior sampling iterations for Stan's Hamiltonian Monte Carlo (HMC) sampler, including warmup.
#' @param adapt_delta A parameter that adjusts the behavior of Stan's Hamiltonian Monte Carlo (HMC) posterior sampling algorithm. Increasing closer to 1 can help avoid divergent transitions.
#' @return A stan model object; check back soon for more information.
#' @details In development; check back soon. See README at github.com/kyle-rosenblad/tradescape.
#' @examples
#'
#' # In development; check back soon.
#'

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
                             cores=cores,
                             chains=chains,
                             warmup=warmup,
                             iter=iter,
                             verbose=F,
                             control=list(adapt_delta=adapt_delta))
  return(stanmod)
}

