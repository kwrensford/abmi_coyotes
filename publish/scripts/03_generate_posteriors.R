#===============================================================================
# GENERATE POSTERIOR PREDICTIONS
#
# Produces:
#
#   abundance ~ climate PC1
#   abundance ~ climate PC2
#   abundance ~ human footprint
#
# using posterior draws from fitted models.
#
#===============================================================================

library(tidyverse)
library(posterior)

dir.create(
  "output/predictions",
  recursive = TRUE,
  showWarnings = FALSE
)

draw_files <- list.files(
  "output/draws",
  pattern = "draws.rds$",
  full.names = TRUE
)

for(f in draw_files){
  
  draws <- readRDS(f)
  
  pars <- draws %>%
    select(
      alpha_lambda,
      beta_clim1,
      beta_clim2,
      beta_human,
      beta_interact1,
      beta_interact2
    )
  
  pc1_grid <- seq(
    -3,
    3,
    length.out = 100
  )
  
  pred_list <- lapply(
    pc1_grid,
    function(x){
      
      exp(
        pars$alpha_lambda +
          pars$beta_clim1 * x
      )
      
    }
  )
  
  pred_mat <- do.call(
    cbind,
    pred_list
  )
  
  pred_summary <- tibble(
    
    PC1 = pc1_grid,
    
    median =
      apply(pred_mat, 2, median),
    
    lower =
      apply(pred_mat, 2,
            quantile,
            probs = 0.025),
    
    upper =
      apply(pred_mat, 2,
            quantile,
            probs = 0.975)
    
  )
  
  outfile <- gsub(
    "_draws.rds",
    "_pc1_prediction.rds",
    basename(f)
  )
  
  saveRDS(
    pred_summary,
    file.path(
      "output/predictions",
      outfile
    )
  )
  
}