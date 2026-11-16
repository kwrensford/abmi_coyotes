#Calculating overlap

library(tidyverse)
library(sf)
library(cmdstanr)
library(posterior)

# function to reconstruct lambda_site from posterior draws
predict_lambda <- function(draws, stan_data) {
  
  exp(
    outer(draws$alpha_lambda,
          rep(1, stan_data$I)) +
      
      outer(draws$beta_clim1,
            stan_data$clim_pc1) +
      
      outer(draws$beta_clim2,
            stan_data$clim_pc2) +
      
      outer(draws$beta_human,
            stan_data$human_fp) +
      
      outer(draws$beta_interact1,
            stan_data$clim_pc1 * stan_data$human_fp) +
      
      outer(draws$beta_interact2,
            stan_data$clim_pc2 * stan_data$human_fp)
  )
}

#List model names
model_files <- list.files(
  "output/fits",
  pattern = "_model\\.rds$",
  full.names = TRUE
)

models <- setNames(
  lapply(model_files, readRDS),
  basename(model_files) |>
    str_remove("_model\\.rds$")
)

#Ientify draw files
draw_files <- list.files(
  "output/draws",
  pattern = "_draws\\.rds$",
  full.names = TRUE
)

#Load in draw files
draws <- draw_files %>%
  setNames(
    basename(draw_files) %>%
      str_remove("_draws\\.rds$")
  ) |>
  lapply(readRDS)

#Generate lambdas
lambda_list <- lapply(
  names(draws),
  function(nm) {
    
    predict_lambda(
      draws = draws[[nm]],
      stan_data = models[[nm]]$stan_data
    )
    
  }
)

names(lambda_list) <- names(draws)

#Joint abundance prob
occ_list <- lapply(
  lambda_list[c(
    "coyote_boreal",
    "lynx_boreal",
    "redfox_boreal"
  )],
  function(x) 1 - exp(-x)
)

joint_occurrence <- Reduce(
  `*`,
  occ_list
)

#Summarize posterior estimates
joint_mean <- colMeans(joint_occurrence)

joint_lcl <- apply(
  joint_occurrence,
  2,
  quantile,
  probs = 0.025
)

joint_ucl <- apply(
  joint_occurrence,
  2,
  quantile,
  probs = 0.975
)

#Region Overlap
boreal_species <- c(
  "coyote_boreal",
  "lynx_boreal",
  "redfox_boreal"
)

boreal_occ <- lapply(
  lambda_list[boreal_species],
  function(x) 1 - exp(-x)
)

boreal_joint <- Reduce(`*`, boreal_occ)

covars_boreal$overlap_mean <- colMeans(boreal_joint)

covars_boreal$overlap_lcl <- apply(
  boreal_joint,
  2,
  quantile,
  0.025
)

covars_boreal$overlap_ucl <- apply(
  boreal_joint,
  2,
  quantile,
  0.975
)

#Load in Covariates
covars <- read.csv("data/covars.csv")
covars_boreal <- read.csv("data/covars_boreal.csv")
covars_park <- read.csv("data/covars_park.csv")

boreal_sites <- covars_boreal %>%
  group_by(site) %>%
  slice(1) %>%
  ungroup()

#Append overlap metrics
boreal_sites$overlap_mean <- joint_mean
boreal_sites$overlap_lcl  <- joint_lcl
boreal_sites$overlap_ucl  <- joint_ucl

#Pairwise Abundance Overlap

#Coyote Lynx
coyote_occ <- 1 - exp(-lambda_list$coyote_boreal)
lynx_occ   <- 1 - exp(-lambda_list$lynx_boreal)

coyote_lynx <- coyote_occ * lynx_occ

coyote_lynx_mean <- colMeans(coyote_lynx)

boreal_sites$coyote_lynx <- coyote_lynx_mean

coyote_lynx_abund <- pmin(
  lambda_list$coyote_boreal,
  lambda_list$lynx_boreal
)

boreal_sites$coyote_lynx_abund <-
  colMeans(coyote_lynx_abund)

#Coyote Red Fox
fox_occ <- 1 - exp(-lambda_list$redfox_boreal)

coyote_fox <- coyote_occ * fox_occ

coyote_fox_mean <- colMeans(coyote_fox)

boreal_sites$coyote_fox  <- coyote_fox_mean

coyote_fox_abund <- pmin(
  lambda_list$coyote_boreal,
  lambda_list$redfox_boreal
)

boreal_sites$coyote_fox_abund <-
  colMeans(coyote_fox_abund)

#Red Fox Lynx
lynx_fox <- lynx_occ * fox_occ

lynx_fox_mean <- colMeans(lynx_fox)

boreal_sites$lynx_fox <- lynx_fox_mean

lynx_fox_abund <- pmin(
  lambda_list$lynx_boreal,
  lambda_list$redfox_boreal
)

boreal_sites$lynx_fox_abund <-
  colMeans(lynx_fox_abund)

boreal_sites$three_species <- boreal_sites$overlap_mean

#Map Overlap
library(sf)
library(ggplot2)
library(patchwork)
library(viridis)

# convert to sf
boreal_sf <- st_as_sf(
  boreal_sites,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

# common scale across all panels
zmax <- max(
  boreal_sf$coyote_lynx,
  boreal_sf$coyote_fox,
  boreal_sf$lynx_fox,
  boreal_sf$three_species,
  na.rm = TRUE
)

#Basemap
library(rnaturalearth)
library(rnaturalearthhires)

ab <- neab <- neab <- ne_states(
  country = "canada",
  returnclass = "sf"
)

ab <- subset(ab, name == "Alberta")

make_panel <- function(data, var, title) {
  
  ggplot() +
    geom_sf(
      data = ab,
      fill = "grey95",
      color = "grey60"
    ) +
    geom_sf(
      data = data,
      aes(color = .data[[var]]),
      size = 2.5
    ) +
    scale_color_viridis_c(
      option = "magma",
      limits = c(0, zmax),
      name = "Overlap"
    ) +
    labs(title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "right"
    )
}

p1 <- make_panel(
  boreal_sf,
  "coyote_lynx",
  "A. Coyote–Lynx"
)

p2 <- make_panel(
  boreal_sf,
  "coyote_fox",
  "B. Coyote–Red Fox"
)

p3 <- make_panel(
  boreal_sf,
  "lynx_fox",
  "C. Lynx–Red Fox"
)

p4 <- make_panel(
  boreal_sf,
  "three_species",
  "D. Three-species overlap"
)

fig <-
  (p1 + p2) /
  (p3 + p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig

ggsave(
  "boreal_overlap_4panel.png",
  fig,
  width = 12,
  height = 10,
  dpi = 600
)
