#===============================================================================
# FIT SPECIES-SPECIFIC N-MIXTURE MODELS
#
# Species:
#   Coyote
#   Red fox
#   Canada lynx
#
# Regions:
#   Boreal
#   Parkland
#
#===============================================================================

library(cmdstanr)
library(tidyverse)

source("R/01_nmix_functions.R")

dir.create("output/fits",
           recursive = TRUE,
           showWarnings = FALSE)

dir.create("output/draws",
           recursive = TRUE,
           showWarnings = FALSE)

dir.create("output/summaries",
           recursive = TRUE,
           showWarnings = FALSE)

#------------------------------------------------------------------------------
# Model registry
#------------------------------------------------------------------------------
#
# Every model fitted in the project is represented here.
#
# Adding future species only requires adding another row.
#
#------------------------------------------------------------------------------

model_sets <- tribble(
  
  ~species, ~region, ~count_file, ~covar_file,
  
  "coyote", "boreal",
  "data/coyotecounts_boreal.csv",
  "data/covars_boreal.csv",
  
  "coyote", "park",
  "data/coyotecounts_park.csv",
  "data/covars_park.csv",
  
  "redfox", "boreal",
  "data/redfoxcounts_boreal.csv",
  "data/covars_boreal.csv",
  
  "redfox", "park",
  "data/redfoxcounts_park.csv",
  "data/covars_park.csv",
  
  "lynx", "boreal",
  "data/lynxcounts_boreal.csv",
  "data/covars_boreal.csv"
  
)

#------------------------------------------------------------------------------
# Compile Stan model once
#------------------------------------------------------------------------------

mod <- cmdstan_model(
  "stan/nmix_nb.stan"
)

#------------------------------------------------------------------------------
# Fit all models
#------------------------------------------------------------------------------

for(i in seq_len(nrow(model_sets))){
  
  species <- model_sets$species[i]
  
  region <- model_sets$region[i]
  
  message(
    "Running ",
    species,
    " (",
    region,
    ")"
  )
  
  counts <- as.matrix(
    read.csv(
      model_sets$count_file[i]
    )
  )
  
  covars <- read.csv(
    model_sets$covar_file[i]
  ) %>%
    scale_covariates()
  
  stan_data <- build_stan_data(
    abund = counts,
    covars = covars
  )
  
  fit <- fit_nmix_model(
    stan_data = stan_data,
    model = mod
  )
  
  model_output <- list(
    
    species = species,
    
    region = region,
    
    stan_data = stan_data,
    
    fit = fit,
    
    fit_date = Sys.Date()
    
  )
  
  saveRDS(
    model_output,
    sprintf(
      "output/fits/%s_%s_model.rds",
      species,
      region
    )
  )
  
  posterior_draws <- fit$draws(
    variables = c(
      "alpha_lambda",
      "beta_clim1",
      "beta_clim2",
      "beta_human",
      "beta_interact1",
      "beta_interact2",
      "alpha_p",
      "beta_bait",
      "sigma_cam",
      "phi"
    ),
    format = "draws_df"
  )
  
  saveRDS(
    posterior_draws,
    sprintf(
      "output/draws/%s_%s_draws.rds",
      species,
      region
    )
  )
  
  write.csv(
    fit$summary(),
    sprintf(
      "output/summaries/%s_%s_summary.csv",
      species,
      region
    ),
    row.names = FALSE
  )
  
}