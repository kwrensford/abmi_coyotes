library(vroom)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(forcats)
library(lubridate)
library(MASS)
library(hms)
library(unmarked)
library(camtrapR)
library(overlap)
library(stars)
library(Rcpp)
library(stringr)
library(mgcv)
library(cowplot)
library(gridExtra)
library(ubms)
library(MCMCvis)
library(sjPlot)
library(sjmisc)
library(ggeffects)
library(cmdstanr)
library(rjags)
library(jagsUI)


coyotecounts <- as.matrix(read.csv("data/coyotecounts.csv"))
coyotecounts_boreal <- as.matrix(read.csv("data/coyotecounts_boreal.csv"))
coyotecounts_park <- as.matrix(read.csv("data/coyotecounts_park.csv"))


redfoxcounts <- as.matrix(read.csv("data/redfoxcounts.csv"))
redfoxcounts_boreal <- as.matrix(read.csv("data/redfoxcounts_boreal.csv"))
redfoxcounts_park <- as.matrix(read.csv("data/redfoxcounts_park.csv"))

lynxcounts <- as.matrix(read.csv("data/lynxcounts.csv"))
lynxcounts_boreal <- as.matrix(read.csv("data/lynxcounts_boreal.csv"))

martencounts <- as.matrix(read.csv("data/martencounts.csv"))
martencounts_boreal <- as.matrix(read.csv("data/martencounts.csv"))

fishercounts <- as.matrix(read.csv("data/fishercounts.csv"))
fishercounts_boreal <- as.matrix(read.csv("data/fishercounts_boreal.csv"))

#Load in covariate data
covars <- read.csv("data/covars.csv")
covars_boreal <- read.csv("data/covars_boreal.csv")
covars_park <- read.csv("data/covars_park.csv")

#Scale data for all covariates
covars$PC1<-scale(covars$PC1)
covars$PC2<-scale(covars$PC2)
covars$hf<-scale(covars$hf)
covars$baited <- as.numeric(covars$baited)

covars_boreal$PC1<-scale(covars_boreal$PC1)
covars_boreal$PC2<-scale(covars_boreal$PC2)
covars_boreal$hf<-scale(covars_boreal$hf)
covars_boreal$baited <- as.numeric(covars_boreal$baited)

covars_park$PC1<-scale(covars_park$PC1)
covars_park$PC2<-scale(covars_park$PC2)
covars_park$hf<-scale(covars_park$hf)
covars_park$baited <- as.numeric(covars_park$baited)

#Format data
abund <- as.matrix(coyotecounts)   # cameras × occasions
K      <- ncol(abund)              # number of occasions

site_id   <- covars$site           # site per camera
sites     <- sort(unique(site_id))
I         <- length(sites)         # number of sites

# cameras per site
J_tab <- table(site_id)
J     <- as.integer(J_tab[match(sites, names(J_tab))])
max_J <- max(J)
max_K <- K

y      <- array(0L, dim = c(I, max_J, max_K))
use_y  <- array(0L, dim = c(I, max_J, max_K))
baited_arr <- matrix(0L, nrow = I, ncol = max_J)

clim_pc1 <- numeric(I)
clim_pc2 <- numeric(I)
human_fp <- numeric(I)

for (i in seq_along(sites)) {
  this_site <- sites[i]
  cams_i    <- which(site_id == this_site)   # row indices in covars/coyotecounts
  J_i       <- length(cams_i)
  
  # site-level covariates: same for all cameras at this site
  clim_pc1[i] <- covars$PC1[cams_i[1]]
  clim_pc2[i] <- covars$PC2[cams_i[1]]
  human_fp[i] <- covars$hf[cams_i[1]]
  
  for (j in seq_len(J_i)) {
    row_idx <- cams_i[j]
    
    # camera-level detection covariate
    baited_arr[i, j] <- as.integer(covars$baited[row_idx])
    
    for (k in seq_len(max_K)) {
      val <- abund[row_idx, k]
      if (!is.na(val)) {
        y[i, j, k]     <- as.integer(val)
        use_y[i, j, k] <- 1L
      }
    }
  }
}

#Compile stan data into list
stan_data <- list(
  I       = I,
  J       = J,
  max_J   = max_J,
  max_K   = max_K,
  y       = y,
  use_y   = use_y,
  clim_pc1 = clim_pc1,
  clim_pc2 = clim_pc2,
  human_fp = human_fp,
  baited   = baited_arr,
  K        = 100  # or some upper bound for N
)

mod <- cmdstan_model(
  "C:/Users/Gaynor Lab/Documents/abmi_coyotes/coyote_nmix.stan"
)

str(stan_data, max.level = 1)

stan_data$J

i <- 1
j <- 1
cbind(
  y      = stan_data$y[i, j, ],
  use_y  = stan_data$use_y[i, j, ]
)

head(stan_data$clim_pc1)
head(stan_data$clim_pc2)
head(stan_data$human_fp)

stan_data$baited[1, 1:stan_data$J[2]]

stan_data$K
