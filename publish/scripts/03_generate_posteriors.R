#===============================================================================
# GENERATE POSTERIOR PREDICTIONS
#
# Produces:
#
#   abundance ~ climate PC1:hf
#   abundance ~ climate PC2:hf
#
# using posterior draws from fitted models.
#
#===============================================================================

library(tidyverse)
library(posterior)

generate_interaction_predictions <- function(
    draws,
    covars,
    focal = c("PC1", "PC2"),
    moderator = "hf",
    mod_probs = c(0.25, 0.75),
    n = 100){
  
  focal <- match.arg(focal)
  
  hf_vals <- quantile(
    covars[[moderator]],
    probs = mod_probs,
    na.rm = TRUE
  )
  
  focal_seq <- seq(
    min(covars[[focal]], na.rm = TRUE),
    max(covars[[focal]], na.rm = TRUE),
    length.out = n
  )
  
  build_curve <- function(hf_value,
                          label){
    
    pred_mat <- sapply(
      
      focal_seq,
      
      function(x){
        
        if(focal == "PC1"){
          
          exp(
            draws$alpha_lambda +
              draws$beta_clim1 * x +
              draws$beta_human * hf_value +
              draws$beta_interact1 * x * hf_value
          )
          
        } else {
          
          exp(
            draws$alpha_lambda +
              draws$beta_clim2 * x +
              draws$beta_human * hf_value +
              draws$beta_interact2 * x * hf_value
          )
          
        }
        
      }
      
    )
    
    tibble(
      
      focal = focal_seq,
      
      pred = apply(
        pred_mat,
        2,
        median
      ),
      
      lower = apply(
        pred_mat,
        2,
        quantile,
        probs = 0.025
      ),
      
      upper = apply(
        pred_mat,
        2,
        quantile,
        probs = 0.975
      ),
      
      hf_group = label
      
    )
    
  }
  
  bind_rows(
    
    build_curve(
      hf_vals[1],
      "Low HF"
    ),
    
    build_curve(
      hf_vals[2],
      "High HF"
    )
    
  )
  
}


covars_boreal <- read.csv(
  "data/covars_boreal.csv"
)

covars_park <- read.csv(
  "data/covars_park.csv"
)

covars_boreal$PC1 <- scale(covars_boreal$PC1)
covars_boreal$PC2 <- scale(covars_boreal$PC2)
covars_boreal$hf  <- scale(covars_boreal$hf)

covars_park$PC1 <- scale(covars_park$PC1)
covars_park$PC2 <- scale(covars_park$PC2)
covars_park$hf  <- scale(covars_park$hf)

draw_files <- list.files(
  "output/draws",
  full.names = TRUE,
  pattern = "_draws.rds$"
)

for(f in draw_files){
  
  draws <- readRDS(f)
  
  if(grepl("boreal", f)){
    
    covars <- covars_boreal
    
  } else {
    
    covars <- covars_park
    
  }
  
  pc1_preds <- generate_interaction_predictions(
    draws,
    covars,
    focal = "PC1"
  )
  
  pc2_preds <- generate_interaction_predictions(
    draws,
    covars,
    focal = "PC2"
  )
  
  saveRDS(
    pc1_preds,
    
    file.path(
      "output/predictions",
      
      gsub(
        "_draws.rds",
        "_pc1.rds",
        
        basename(f)
      )
    )
  )
  
  saveRDS(
    pc2_preds,
    
    file.path(
      "output/predictions",
      
      gsub(
        "_draws.rds",
        "_pc2.rds",
        
        basename(f)
      )
    )
  )
  
}

}