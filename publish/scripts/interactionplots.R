#Building Interaction Plots

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
library(ggimage)
library(cowplot)
library(patchwork)
library(grid)
library(magick)
library(tibble)


#Load in Model Outputs

coyote_nmix_ubms_boreal <- readRDS("coyote_nmix_ubms_boreal.rds")
coyote_nmix_ubms_park <- readRDS("coyote_nmix_ubms_park.rds")
redfox_nmix_stan_boreal <- readRDS("redfox_nmix_stan_boreal.rds")
redfox_nmix_stan_park <- readRDS("redfox_nmix_stan_park.rds")
lynx_nmix_stan <- readRDS("lynx_nmix_stan.rds")

#Load in Covariates
covars <- read.csv("data/covars.csv")
covars_boreal <- read.csv("data/covars_boreal.csv")
covars_park <- read.csv("data/covars_park.csv")

###Predict interaction plot function
plot_interaction_ubms <- function(model,
                                  focal,
                                  moderator,
                                  species = NULL,
                                  submodel = "state",
                                  data,
                                  mod_probs = c(0.25, 0.75),
                                  n = 1000,
                                  ylab = NULL,
                                  xlab = NULL,
                                  title = NULL,
                                  region = NULL,
                                  ylim = NULL) {
  
  library(ggplot2)
  library(colorspace)

  
  # Region color base
  region_colors <- list(
    boreal = "#1b9e77",
    park   = "#d95f02"
  )
  base_col <- region_colors[[region]]
  col_low  <- lighten(base_col, 0.4)
  col_high <- darken(base_col, 0.4)
  
  # Moderator quantiles
  mod_vals <- quantile(data[[moderator]], probs = mod_probs, na.rm = TRUE)
  
  # Build prediction grids
  make_grid <- function(mod_value) {
    grid <- data.frame(
      focal = seq(min(data[[focal]], na.rm = TRUE),
                  max(data[[focal]], na.rm = TRUE),
                  length.out = n)
    )
    names(grid)[1] <- focal
    
    grid[[moderator]] <- mod_value
    
    numeric_covs <- names(data)[sapply(data, is.numeric)]
    other_covs <- setdiff(numeric_covs, c(focal, moderator))
    for (v in other_covs) {
      grid[[v]] <- mean(data[[v]], na.rm = TRUE)
    }
    grid
  }
  
  grid_low  <- make_grid(mod_vals[1])
  grid_high <- make_grid(mod_vals[2])
  
  # Predictions
  pred_low <- ubms::predict(model, submodel = submodel,
                            species = species, newdata = grid_low)
  pred_high <- ubms::predict(model, submodel = submodel,
                             species = species, newdata = grid_high)
  
  # Long-format dataframe for legend
  df_long <- rbind(
    data.frame(
      cov = grid_low[[focal]],
      pred = pred_low$Predicted,
      lower = pred_low$`2.5%`,
      upper = pred_low$`97.5%`,
      hf_group = "Low HF"
    ),
    data.frame(
      cov = grid_high[[focal]],
      pred = pred_high$Predicted,
      lower = pred_high$`2.5%`,
      upper = pred_high$`97.5%`,
      hf_group = "High HF"
    )
  )
  
  # Labels
  if (is.null(ylab)) ylab <- paste("Predicted", species)
  if (is.null(xlab)) xlab <- NULL  # we annotate manually
  
  # Plot
  p <- ggplot(df_long, aes(x = cov, y = pred,
                           color = hf_group, fill = hf_group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.25, color = NA) +
    geom_line(size = 1) +
    
    scale_color_manual(values = c("Low HF" = col_low, "High HF" = col_high),
                       name = "Human Footprint") +
    scale_fill_manual(values = c("Low HF" = col_low, "High HF" = col_high),
                      name = "Human Footprint") +
    labs(y = ylab, x = NULL, title = title) +
    theme_classic() +
    theme(legend.position = "none") 
  
  # Apply y-limits if provided
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}


#Calculate y limits for each species/region to ensure consistent y axis
get_ylim <- function(model, focal, moderator, species, data, mod_probs = c(0.25, 0.75), n = 1000) {
  
  mod_vals <- quantile(data[[moderator]], probs = mod_probs, na.rm = TRUE)
  
  make_grid <- function(mod_value) {
    grid <- data.frame(
      focal = seq(min(data[[focal]], na.rm = TRUE),
                  max(data[[focal]], na.rm = TRUE),
                  length.out = n)
    )
    names(grid)[1] <- focal
    
    grid[[moderator]] <- mod_value
    
    numeric_covs <- names(data)[sapply(data, is.numeric)]
    other_covs <- setdiff(numeric_covs, c(focal, moderator))
    for (v in other_covs) {
      grid[[v]] <- mean(data[[v]], na.rm = TRUE)
    }
    grid
  }
  
  grid_low  <- make_grid(mod_vals[1])
  grid_high <- make_grid(mod_vals[2])
  
  pred_low  <- ubms::predict(model, submodel = "state", species = species, newdata = grid_low)
  pred_high <- ubms::predict(model, submodel = "state", species = species, newdata = grid_high)
  
  range(c(pred_low$`2.5%`, pred_low$`97.5%`,
          pred_high$`2.5%`, pred_high$`97.5%`),
        na.rm = TRUE)
}

ylim_coyote_boreal <- get_ylim(coyote_nmix_ubms_boreal, "PC1", "hf", "coyote", covars_boreal)
ylim_coyote_park <- get_ylim(coyote_nmix_ubms_park, "PC1", "hf", "coyote", covars_boreal)
ylim_redfox_boreal <- get_ylim(redfox_nmix_stan_boreal, "PC1", "hf", "redfox", covars_boreal)
ylim_redfox_park <- get_ylim(redfox_nmix_stan_park, "PC1", "hf", "redfox", covars_boreal)
ylim_lynx_boreal <- get_ylim(lynx_nmix_stan, "PC1", "hf", "lynx", covars_boreal)


#Create Plots for all species/region models
plots <- list(
  coyote_boreal_pc1 = plot_interaction_ubms(
    model = coyote_nmix_ubms_boreal,
    focal = "PC1", moderator = "hf",
    species = "coyote", data = covars_boreal,
    region = "boreal",
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL, ylim = ylim_coyote_boreal
  ),
  
  coyote_boreal_pc2 = plot_interaction_ubms(
    model = coyote_nmix_ubms_boreal,
    focal = "PC2", moderator = "hf",
    species = "coyote", data = covars_boreal,
    region = "boreal",
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL, ylim = ylim_coyote_boreal 
  ),
  coyote_park_pc1 = plot_interaction_ubms(
    model = coyote_nmix_ubms_park,
    focal = "PC1", moderator = "hf",
    species = "coyote", data = covars_park,
    region = "park",
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL, ylim = ylim_coyote_park
  ),
  
  coyote_park_pc2 = plot_interaction_ubms(
    model = coyote_nmix_ubms_park,
    focal = "PC2", moderator = "hf",
    species = "coyote", data = covars_park,
    region = "park",
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL, ylim = ylim_coyote_park
  ),
  
  redfox_boreal_pc1 = plot_interaction_ubms(
    model = redfox_nmix_stan_boreal,
    focal = "PC1", moderator = "hf",
    species = "redfox", data = covars_boreal,
    region = "boreal",
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL, ylim = ylim_redfox_boreal
  ),
  
  redfox_boreal_pc2 = plot_interaction_ubms(
    model = redfox_nmix_stan_boreal,
    focal = "PC2", moderator = "hf",
    species = "redfox", data = covars_boreal,
    region = "boreal",
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL, ylim = ylim_redfox_boreal
  ),
  
  redfox_park_pc1 = plot_interaction_ubms(
    model = redfox_nmix_stan_park,
    focal = "PC1", moderator = "hf",
    species = "redfox", data = covars_park,
    region = "park",
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL, ylim  = ylim_redfox_park
  ),
  
  redfox_park_pc2 = plot_interaction_ubms(
    model = redfox_nmix_stan_park,
    focal = "PC2", moderator = "hf",
    species = "redfox", data = covars_park,
    region = "park",
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL, ylim = ylim_redfox_park
  ),
  
  lynx_boreal_pc1 = plot_interaction_ubms(
    model = lynx_nmix_stan,
    focal = "PC1", moderator = "hf",
    species = "lynx", data = covars_boreal,
    region = "boreal",
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL, ylim = ylim_lynx_boreal
  ),
  
  lynx_boreal_pc2 = plot_interaction_ubms(
    model = lynx_nmix_stan,
    focal = "PC2", moderator = "hf",
    species = "lynx", data = covars_boreal,
    region = "boreal",
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL, ylim = ylim_lynx_boreal
  )
)

#Function to add species icon
make_icon_grob <- function(species_name, size = .75) {
  icon_path <- file.path("data/phylopics", paste0(species_name, ".png"))
  img <- magick::image_read(icon_path)
  rasterGrob(img, width = size, height = size)
}


#Create x axis labels
xlab_pc1 <- ggplot() +
  annotate("text", x = 0.5, y = 0.65,
           label = "Cold → Warm gradient",
           size = 3.2, color = "grey40", lineheight = 0.9) +
  annotate("text", x = 0.5, y = 0.35,
           label = "Climate PC1",
           size = 5, fontface = "bold", lineheight = 0.9) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))


xlab_pc2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.65,
           label = "Wet → Dry gradient",
           size = 3.2, color = "grey40", lineheight = 0.9) +
  annotate("text", x = 0.5, y = 0.35,
           label = "Climate PC2",
           size = 5, fontface = "bold", lineheight = 0.9) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))




#Create a left side label strip for each row, comine icon and label vertically

make_row_label <- function(species_name, region_name) {
  
  icon_path <- file.path("data/phylopics", paste0(species_name, ".png"))
  
  icon_panel <- ggplot() +
    geom_image(aes(x = 0.5, y = 0.5, image = icon_path),
               size = 1.75,
               asp = TRUE) +
    theme(
      plot.margin = margin(4, 0, 4, 0)
    ) + 
    theme_void()
  
  species_panel <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = species_name,
             size = 5, fontface = "bold") +
    theme_void()
  
  region_panel <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = region_name,
             size = 5) +
    theme_void()
  
  # Combine vertically
  combined <- icon_panel / species_panel / region_panel
  
  # rap in fixed-height panel so patchwork doesn't stretch the row
  wrap_elements(full = combined, clip = TRUE)
}

#Combine row label + PC1 + PC2 into a row
row_coyote_park <- 
  make_row_label("Coyote", "Grassland/Parkland") |
  plots$coyote_park_pc1 |
  plots$coyote_park_pc2

row_coyote_boreal <- 
  make_row_label("Coyote", "Boreal/Shield") |
  plots$coyote_boreal_pc1 |
  plots$coyote_boreal_pc2

row_redfox_park <- 
  make_row_label("Red fox", "Grassland/Parkland") | 
  plots$redfox_park_pc1 | 
  plots$redfox_park_pc2

row_redfox_boreal <- 
  make_row_label("Red fox", "Boreal/Shield") | 
  plots$redfox_boreal_pc1 | 
  plots$redfox_boreal_pc2

row_lynx_boreal <- 
  make_row_label("Canada lynx",   "Boreal/Shield") | 
  plots$lynx_boreal_pc1 | 
  plots$lynx_boreal_pc2

#Stack rows into the final grid
grid <- 
  row_coyote_park /
  row_coyote_boreal /
  row_redfox_park /
  row_redfox_boreal /
  row_lynx_boreal

# Apply column widths BEFORE adding PC labels
# left strip = narrow, PC1 = wide, PC2 = wide
grid <- grid + plot_layout(widths = c(1, 3, 3))

#Build region gradient bar legend

make_hf_bar <- function(region_name, low_col, high_col) {
  
  grad_df <- tibble(
    x = seq(0, 1, length.out = 200),
    y = 1,
    col = colorRampPalette(c(low_col, high_col))(200)
  )
  
  ggplot(grad_df, aes(x, y, fill = col)) +
    geom_tile(height = 0.2) +
    scale_fill_identity() +
    
    annotate("text", x = 0.5, y = 1.3,
             label = region_name,
             size = 4, fontface = "bold") +
    
    annotate("text", x = 0, y = 0.4,
             label = "Low HF", size = 3.2, hjust = 0) +
    annotate("text", x = 1, y = 0.4,
             label = "High HF", size = 3.2, hjust = 1) +
    
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(2, 2, 2, 2)
    )
}

boreal_bar <- make_hf_bar(
  "Boreal/Shield",
  lighten(region_colors$boreal, 0.4),
  darken(region_colors$boreal, 0.4)
)

park_bar <- make_hf_bar(
  "Grassland/Parkland",
  lighten(region_colors$park, 0.4),
  darken(region_colors$park, 0.4)
)


hf_region_legend <- boreal_bar / park_bar

pc_labels_row <- hf_region_legend | xlab_pc1 | xlab_pc2

grid <- grid / pc_labels_row


#Combine grid + legend with controlled widths
grid <- grid + plot_layout(widths = c(1,3,3), 
                                             heights = c(1,1,1,1, 1, .25))

final <- grid


final

ggsave("interactionplots_abmi.pdf", final, width = 7, height = 12)
ggsave("interactionplots_abmi.png", final, width = 7, height = 12)
