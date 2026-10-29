#===============================================================================
# FUNCTIONS FOR NEGATIVE BINOMIAL N-MIXTURE ANALYSES
#
# Purpose:
#   Utility functions used across model-fitting and prediction workflows.
#
# Author:
#  Kwasi Wrensford
#
#===============================================================================

library(dplyr)

#------------------------------------------------------------------------------
# Scale covariates used in abundance and detection models
#
# Args:
#   covars : data frame containing site-level and camera-level covariates
#
# Returns:
#   Data frame with standardized predictors
#------------------------------------------------------------------------------

scale_covariates <- function(covars){
  
  covars %>%
    mutate(
      
      PC1 = as.numeric(scale(PC1)),
      
      PC2 = as.numeric(scale(PC2)),
      
      hf = as.numeric(scale(hf)),
      
      baited = as.integer(baited)
      
    )
  
}

#------------------------------------------------------------------------------
# Build Stan-formatted data object
#
# Converts camera-level count matrices into a site × camera × occasion array.
#
# Dimensions:
#
#   i = site
#   j = camera
#   k = sampling occasion
#
# Returns:
#   Named list suitable for cmdstanr
#------------------------------------------------------------------------------

build_stan_data <- function(abund,
                            covars,
                            abundance_upper = 100){
  
  site_id <- covars$site
  
  sites <- sort(unique(site_id))
  
  I <- length(sites)
  
  K_occasions <- ncol(abund)
  
  J_tab <- table(site_id)
  
  J <- as.integer(
    J_tab[match(sites, names(J_tab))]
  )
  
  max_J <- max(J)
  
  y <- array(
    0L,
    dim = c(I, max_J, K_occasions)
  )
  
  use_y <- array(
    0L,
    dim = c(I, max_J, K_occasions)
  )
  
  baited_arr <- matrix(
    0L,
    nrow = I,
    ncol = max_J
  )
  
  clim_pc1 <- numeric(I)
  clim_pc2 <- numeric(I)
  human_fp <- numeric(I)
  
  for(i in seq_along(sites)){
    
    cams_i <- which(site_id == sites[i])
    
    clim_pc1[i] <- covars$PC1[cams_i[1]]
    clim_pc2[i] <- covars$PC2[cams_i[1]]
    human_fp[i] <- covars$hf[cams_i[1]]
    
    for(j in seq_along(cams_i)){
      
      row_idx <- cams_i[j]
      
      baited_arr[i,j] <- as.integer(
        covars$baited[row_idx]
      )
      
      for(k in seq_len(K_occasions)){
        
        val <- abund[row_idx,k]
        
        if(!is.na(val)){
          
          y[i,j,k] <- as.integer(val)
          
          use_y[i,j,k] <- 1L
          
        }
      }
    }
  }
  
  list(
    
    I = I,
    
    J = J,
    
    max_J = max_J,
    
    max_K = K_occasions,
    
    y = y,
    
    use_y = use_y,
    
    clim_pc1 = clim_pc1,
    
    clim_pc2 = clim_pc2,
    
    human_fp = human_fp,
    
    baited = baited_arr,
    
    K = abundance_upper
    
  )
  
}

#------------------------------------------------------------------------------
# Fit Stan N-mixture model
#
# Args:
#   stan_data : output from build_stan_data()
#   model     : compiled cmdstanr model
#
# Returns:
#   CmdStanMCMC object
#------------------------------------------------------------------------------

fit_nmix_model <- function(stan_data,
                           model){
  
  model$sample(
    
    data = stan_data,
    
    chains = 4,
    
    parallel_chains = 4,
    
    iter_warmup = 2000,
    
    iter_sampling = 20000,
    
    adapt_delta = 0.95,
    
    max_treedepth = 12,
    
    refresh = 100
    
  )
  
}