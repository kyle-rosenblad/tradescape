// exponentiated quadratic covariance decay function for spatial Gaussian process
functions {
  matrix gp_exp_quad(matrix x, real etasq, real rhosq, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    for (i in 1:N) {
      K[i, i] = etasq + delta;
      for (j in (i+1):N) {
        real sq_dist = 0;
        for (d in 1:2) {
          sq_dist += ((x[i,d] - x[j,d])^2);
        }
        K[i, j] = etasq * exp(-rhosq*sq_dist);
        K[j, i] = K[i, j];
      }
    }
    return K;
  }
}

data {
  // dimensions
  int<lower=0> N_snps;  // number of SNPs
  int<lower=0> N_obs; // total number of observations
  int<lower=0> N_sites; // number of unique sites
  int<lower=0> N_envars; // number of environmental variables
  int<lower=0> ploidy; // maximum allele count for each observation
  
  // indices
  array[N_obs] int<lower=1, upper=N_snps> snp_idx;    // SNP index for each observation
  array[N_obs] int<lower=1, upper=N_sites> site_idx;  // site index for each observation
  
  // data
  array[N_obs] int<lower=0, upper=ploidy> allele_count;    // response variable, allele frequency
  matrix[N_obs, N_envars] envdata;                           // environment data
  matrix[N_sites, 2] coords;                          // [x,y] spatial coordinates for each site
}

parameters {
  // full dataset-level parameters for non-spatial effects
  real<lower=0> sigma_intercept;          // sd of SNP intercepts
  cholesky_factor_corr[N_envars] L_pop;          // Cholesky for environment effects correlation
  vector<lower=0>[N_envars] sigma_pop;           // sd for environment effects
  
  // spatial Gaussian process parameters shared across SNPs
  real<lower=0> etasq;                      // greater values indicate faster decay of covariance with distance
  real<lower=0> rhosq;                    // greater values indicate covariance decay starts from a greater covariance value for nearby locations
  real<lower=0> delta;                    // additive offset for diagonals of variance-covariance matrix; i.e., "intra-location variance"
  
  // SNP-level parameters
  vector[N_snps] intercept_snp;               // intercepts
  array[N_snps] vector[N_envars] z_environment;      // standardized environment effects
  
  // spatial Gaussian process draws (non-centered parameterization)
  array[N_snps] vector[N_sites] z_gp;
}

transformed parameters {
  // transform environment effects
  corr_matrix[N_envars] R_pop = multiply_lower_tri_self_transpose(L_pop);
  cov_matrix[N_envars] Sigma_pop = quad_form_diag(R_pop, sigma_pop);
  array[N_snps] vector[N_envars] beta_environment;
  
  // store spatial Gaussian process draws in one big array
  array[N_snps] vector[N_sites] f;  // spatial effects
  
  // tranform back to centered parameterization for environment effects
  for (i in 1:N_snps) {
    beta_environment[i] = L_pop * (sigma_pop .* z_environment[i]);
  }
  
  // spatial covariance matrix & cholesky
  matrix[N_sites, N_sites] K = gp_exp_quad(coords, etasq, rhosq, delta);
  matrix[N_sites, N_sites] L_K = cholesky_decompose(K);
  
  for (i in 1:N_snps) {
    f[i] = L_K * z_gp[i];
  }
}

model {
  // full dataset-level priors
  sigma_intercept ~ normal(0,5);
  sigma_pop ~ normal(0,5);
  L_pop ~ lkj_corr_cholesky(1);
  rhosq ~ gamma(5, 2);        // prior for rhosq; greater values mean faster decay of covariance with distance
  etasq ~ gamma(5, 5);    // prior etasq; greater values mean covariance decay starts from a greater covariance value for nearby locations
  delta ~ gamma(5, 5);        // additive offset for diagonals of variance-covariance matrix; i.e., "intra-location variance"
  
  // SNP-level priors
  // posteriors can sometimes come out bimodal due to reference allele coin-flip
  // if major allele counts tend to be very high.
  // testing suggests it's not a problem for inference on parameters of interest,
  // but maybe try a symmetric mixture of two Gaussians in future release to more elegantly account for this.
  intercept_snp ~ normal(0, sigma_intercept); 
  for (i in 1:N_snps) {
    z_environment[i] ~ std_normal();
    z_gp[i] ~ std_normal();
  }
  
  // likelihood
  {
    vector[N_obs] logit_p;
    for (n in 1:N_obs) {
      logit_p[n] = intercept_snp[snp_idx[n]] + 
                   dot_product(beta_environment[snp_idx[n]], envdata[n,]') +
                   f[snp_idx[n], site_idx[n]];
    }
    allele_count ~ binomial_logit(ploidy, logit_p);
  }
}
