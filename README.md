**Hierarchical Bayesian modeling of local adaptation to multidimensional environmental gradients with SNP data.**

To install, first follow the rstan installation instructions: (https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). Then, in R, run install.packages(“devtools”), then devtools::install_github(“kyle-rosenblad/tradescape”).

In development. Please don't hesitate to reach out with ideas or bug reports!

A core question in ecology and evolutionary biology is how multidimensional environmental gradients shape evolutionary adaptation. For example, do genetic variants that increase fitness in dry conditions tend to increase or decrease fitness in cold conditions? The answers hold important implications for conservation strategies like assisted gene flow under climate change. E.g., does importing drought-resistant genotypes increase or decrease population vulnerability to cold snaps?

tradescape is designed for exploring these types of questions with biallelic SNP data--i.e., genomic data for some number of individuals across some number of geographic sampling sites, such that some individuals have one allele at a given point in the genome, and other individuals have a different allele. In this scenario, we can cautiously infer patterns of natural selection through principled modeling of the relationship between environmental conditions and allele frequencies. The core assumption is that, if we've structured our model well, associations between allele frequencies and environmental gradients reflect selection patterns.

**Semi-technical summary**
The response variable is the count of the reference allele for one SNP in one individual. We can think of each SNP as having its own sub-model. This sub-model is essentially a logit GLMM with a binomial response variable, which is the mean frequency of the reference allele for the focal SNP. (A beta-binomial version is coming soon. These are helpful for handling certain technical problems that sometimes arise in binomial models.) The key terms in the linear predictor describe the effects of our environmental variables of interest (e.g., temperature, precipitation, etc.) on the reference allele frequency for the focal SNP. Each environmental variable has a linear effect size (on the logit scale) indicating how strongly the mean reference allele frequency increases or decreases with the that variable. Other terms in the linear predictor include a SNP-level intercept (these are drawn from a shared distribution across SNPs) and a spatial Gaussian process--i.e., a flexibly shaped smooth function--designed to account for what we might call "isolation by distance" or "neutral demographic history".

To answer our core questions--e.g., "Do alleles that increase with temperature tend to increase or decrease with precipitation?"--we have to zoom out to the level of the whole data set. At this level, we are modeling associations among the effect sizes for different environmental variables across SNPs. For example, if alleles that increase in frequency with temperature tend to decrease in frequency with precipitation, then our posterior estimates would show a negative correlation between the SNP-level effect sizes for these two environmental variables. (These assocations are modeled using a multivariate normal distribution, from which we can derive a correlation matrix.)

**Full model structure**
**Data Structure**

**Environment variables**: indexed by *m* = 1, ..., *M*

**SNPs**: indexed by *i* = 1, ..., *I*

**Individuals**: indexed by *j* = 1, ..., *J*

**Sites**: indexed by *k* = 1, ..., *K*

Individuals are nested within sites: *j* ∈ *k*

**Response Variable**

For each SNP *i* and individual *j*:*y*ᵢⱼ ~ Binomial(4, *p*ᵢⱼ)

where *y*ᵢⱼ is the count of reference alleles (out of 4 trials) and
*p*ᵢⱼ is the reference allele frequency.

**Linear Predictor**

The linear predictor for the logit-transformed probability:

logit(*p*ᵢⱼ) = Σₘ₌₁^M^ *β*ᵢₘ·*x*ₘₖ + *α*ᵢ + *γ*ᵢₖ

where:

*x*ₘₖ: the *m*-th environment variable at site *k*

*β*ᵢₘ: the effect of the *m*-th environment variable for SNP *i*

*α*ᵢ: the random intercept for SNP *i*

*γ*ᵢₖ: the spatial random effect for SNP *i* at site *k*

**Random Effects**

**SNP-level random intercepts**: *α*ᵢ ~ N(0, *σ*²), *i* = 1, ..., *I*

**Spatial random effects**: For each SNP *i*, the vector of spatial
effects *γ*ᵢ = (*γ*ᵢ₁, *γ*ᵢ₂, ..., *γ*ᵢₖ)ᵀ follows: *γ*ᵢ ~MVN(0,
**K**) where **K** is the covariance matrix with elements: *K*ₖₖ' =
*η*² · exp(-*ρ*² · *D*²ₖₖ') + *δ* · *I*ₖₖ'*D*ₖₖ': spatial distance
between sites *k* and *k'*

*η*²: maximum covariance (when distance = 0)

*ρ*²: rate of covariance decay with distance

*δ*: additional variance for diagonal elements

*I*ₖₖ': indicator function (1 if *k* = *k'*, 0 otherwise)

**Environment effect coefficients**: For each SNP *i*, the vector of
environment effects *β*ᵢ = (*β*ᵢ₁, *β*ᵢ₂, ..., *β*ᵢ~M~)ᵀ follows:*β*ᵢ
~ MVN(0, **Σ**) where **Σ** is a M × M covariance matrix to be
estimated.

**Parameters to Estimate for Random Effects**

*σ*²: variance of SNP-level random intercepts

*η*²: maximum spatial covariance

*ρ*²: spatial covariance decay rate

*δ*: additive offset in diagonal elements of spatial covariance matrix

**Σ**: M × M covariance matrix for environment effect coefficients
