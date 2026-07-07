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

#Load in Model Outputs

coyote_nmix_ubms_boreal <- readRDS("coyote_nmix_ubms_boreal.rds")
coyote_nmix_ubms_park <- readRDS("coyote_nmix_ubms_park.rds")
redfox_nmix_stan_boreal <- readRDS("redfox_nmix_stan_boreal")
redfox_nmix_stan_park <- readRDS("redfox_nmix_stan_park.rds")
lynx_nmix_stan <- readRDS("lynx_nmix_stan.rds")

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
                                  col_low = "darkgreen",
                                  col_high = "darkorange3") {
  
  # 1. Extract moderator quantiles
  mod_vals <- quantile(data[[moderator]], probs = mod_probs, na.rm = TRUE)
  
  # 2. Build prediction grids
  make_grid <- function(mod_value) {
    grid <- data.frame(
      focal = seq(min(data[[focal]], na.rm = TRUE),
                  max(data[[focal]], na.rm = TRUE),
                  length.out = n)
    )
    names(grid)[1] <- focal
    
    # Add moderator
    grid[[moderator]] <- mod_value
    
    # Fix all other covariates at their means
    numeric_covs <- names(data)[sapply(data, is.numeric)]
    other_covs <- setdiff(numeric_covs, c(focal, moderator))
    for (v in other_covs) {
      grid[[v]] <- mean(data[[v]], na.rm = TRUE)
    }
    grid
  }
  
  grid_low <- make_grid(mod_vals[1])
  grid_high <- make_grid(mod_vals[2])
  
  # 3. Predict
  pred_low <- ubms::predict(model,
                            submodel = submodel,
                            species = species,
                            newdata = grid_low)
  
  pred_high <- ubms::predict(model,
                             submodel = submodel,
                             species = species,
                             newdata = grid_high)
  
  # 4. Combine into a plotting dataframe
  df <- data.frame(
    cov = grid_low[[focal]],
    pred_lo = pred_low$Predicted,
    lower_lo = pred_low$`2.5%`,
    upper_lo = pred_low$`97.5%`,
    pred_hi = pred_high$Predicted,
    lower_hi = pred_high$`2.5%`,
    upper_hi = pred_high$`97.5%`
  )
  
  # 5. Labels
  if (is.null(ylab)) ylab <- paste("Predicted", species)
  if (is.null(xlab)) xlab <- focal
  #if (is.null(title)) title <- paste("Effect of", focal, "x", moderator)
  
  # 6. Plot
  library(ggplot2)
  ggplot(df, aes(x = cov)) +
    geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo),
                fill = col_low, alpha = 0.2) +
    geom_line(aes(y = pred_lo), color = col_low, linewidth = 1) +
    geom_ribbon(aes(ymin = lower_hi, ymax = upper_hi),
                fill = col_high, alpha = 0.2) +
    geom_line(aes(y = pred_hi), color = col_high, linewidth = 1) +
    labs(y = ylab, x = xlab, title = title) +
    theme_classic() + 
    theme(legend.position = "none")
}


#Create Plots for all species/region models
plots <- list(
  coyote_boreal_pc1 = plot_interaction_ubms(
    model = coyote_nmix_ubms_boreal,
    focal = "PC1", moderator = "hf",
    species = "coyote", data = covars_boreal,
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL
  ),
  
  coyote_boreal_pc2 = plot_interaction_ubms(
    model = coyote_nmix_ubms_boreal,
    focal = "PC2", moderator = "hf",
    species = "coyote", data = covars_boreal,
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL
  ),
  coyote_park_pc1 = plot_interaction_ubms(
    model = coyote_nmix_ubms_park,
    focal = "PC1", moderator = "hf",
    species = "coyote", data = covars_park,
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL
  ),
  
  coyote_park_pc2 = plot_interaction_ubms(
    model = coyote_nmix_ubms_park,
    focal = "PC2", moderator = "hf",
    species = "coyote", data = covars_park,
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL
  ),
  
  redfox_boreal_pc1 = plot_interaction_ubms(
    model = redfox_nmix_stan_boreal,
    focal = "PC1", moderator = "hf",
    species = "redfox", data = covars_boreal,
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL
  ),
  
  redfox_boreal_pc2 = plot_interaction_ubms(
    model = redfox_nmix_stan_boreal,
    focal = "PC2", moderator = "hf",
    species = "redfox", data = covars_boreal,
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL
  ),
  
  redfox_park_pc1 = plot_interaction_ubms(
    model = redfox_nmix_stan_park,
    focal = "PC1", moderator = "hf",
    species = "redfox", data = covars_park,
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL
  ),
  
  redfox_park_pc2 = plot_interaction_ubms(
    model = redfox_nmix_stan_park,
    focal = "PC2", moderator = "hf",
    species = "redfox", data = covars_park,
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL
  ),
  
  lynx_boreal_pc1 = plot_interaction_ubms(
    model = lynx_nmix_stan,
    focal = "PC1", moderator = "hf",
    species = "lynx", data = covars_boreal,
    ylab = "Abundance", xlab = "Climate PC1",
    title = NULL
  ),
  
  lynx_boreal_pc2 = plot_interaction_ubms(
    model = lynx_nmix_stan,
    focal = "PC2", moderator = "hf",
    species = "lynx", data = covars_boreal,
    ylab = "Abundance", xlab = "Climate PC2",
    title = NULL
  )
)

#Function to add species icon
make_icon_grob <- function(species_name, size = 0.30) {
  icon_path <- file.path("data/phylopics", paste0(species_name, ".png"))
  img <- magick::image_read(icon_path)
  rasterGrob(img, width = size, height = size)
}

#Make a region label panel
make_region_label <- function(region_name) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = region_name, size = 5) +
    theme_void()
}

#Create a left side label strip for each row, comine icon and label vertically
make_row_label <- function(species_name, region_name) {
  icon_panel <- ggplot() +
    annotation_custom(make_icon_grob(species_name),
                      xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = Inf) +
    theme_void() +
    coord_cartesian(clip = "off")
  
  region_panel <- make_region_label(region_name)
  
  icon_panel / region_panel   # vertical stack
}

#Combine row label + PC1 + PC2 into a row
row_coyote_park <- 
  make_row_label("coyote", "Grassland & Parkland") |
  plots$coyote_park_pc1 |
  plots$coyote_park_pc2

row_coyote_boreal <- 
  make_row_label("coyote", "Boreal & Shield") |
  plots$coyote_boreal_pc1 |
  plots$coyote_boreal_pc2

row_redfox_park <- 
  make_row_label("redfox", "Grassland & Parkland") | 
  plots$redfox_park_pc1 | 
  plots$redfox_park_pc2

row_redfox_boreal <- 
  make_row_label("redfox", "Boreal & Shield") | 
  plots$redfox_boreal_pc1 | 
  plots$redfox_boreal_pc2

row_lynx_boreal <- 
  make_row_label("lynx",   "Boreal & Shield") | 
  plots$lynx_boreal_pc1 | 
  plots$lynx_boreal_pc2

#Stack rows into the final grid
grid <- 
  row_coyote_park /
  row_coyote_boreal /
  row_redfox_park /
  row_redfox_boreal /
  row_lynx_boreal

#Remove x axis labels from all but the bottom row
final <- grid +
  plot_annotation(
    tag_levels = NULL,
    theme = theme(
      plot.title = element_text(size = 14)
    )
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

#Make a PC1 axis label panel
xlab_pc1 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Climate PC1", size = 5) +
  theme_void()

#Make a PC2 axis label panel
xlab_pc2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Climate PC2", size = 5) +
  theme_void()

#Attach them under grid
final <- final /
  (xlab_pc1 | xlab_pc2)



#Add column level x axis labels
final <- grid +
  plot_annotation(
    tag_levels = NULL,
    theme = theme(
      plot.title = element_text(size = 14)
    )
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

final <- wrap_plots(
  grid,
  widths = c(1, 3, 3),
  guides = "collect"
) +
  plot_annotation(
    title = NULL,
    caption = NULL
  )

#Add shared legend for human footprint
legend <- get_legend(
  plots$coyote_park_pc1 + theme(legend.position = "right")
)

final <- grid | legend

final

ggsave("interactionplots_abmi.pdf", final, width = 7, height = 12)
ggsave("interactionplots_abmi.png", final, width = 7, height = 12)