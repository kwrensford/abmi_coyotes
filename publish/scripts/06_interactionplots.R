#Building Interaction Plots

library(vroom)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggeffects)
library(ggimage)
library(cowplot)
library(patchwork)
library(grid)
library(magick)
library(tibble)
library(ggtext)


#Load in Predictions 
coyote_boreal_pc1 <- readRDS(
  "output/predictions/coyote_boreal_pc1.rds"
)

coyote_boreal_pc2 <- readRDS(
  "output/predictions/coyote_boreal_pc2.rds"
)

coyote_park_pc1 <- readRDS(
  "output/predictions/coyote_park_pc1.rds"
)

coyote_park_pc2 <- readRDS(
  "output/predictions/coyote_park_pc2.rds"
)

redfox_boreal_pc1 <- readRDS(
  "output/predictions/redfox_boreal_pc1.rds"
)

redfox_boreal_pc2 <- readRDS(
  "output/predictions/redfox_boreal_pc2.rds"
)

redfox_park_pc1 <- readRDS(
  "output/predictions/redfox_park_pc1.rds"
)

redfox_park_pc2 <- readRDS(
  "output/predictions/redfox_park_pc2.rds"
)

lynx_boreal_pc1 <- readRDS(
  "output/predictions/lynx_boreal_pc1.rds"
)

lynx_boreal_pc2 <- readRDS(
  "output/predictions/lynx_boreal_pc2.rds"
)

#Load in Covariates
covars <- read.csv("data/covars.csv")
covars_boreal <- read.csv("data/covars_boreal.csv")
covars_park <- read.csv("data/covars_park.csv")

#Scale covariates
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

#Calculate global limits for PC1 and PC2
pc1_lim <- range(
  covars_boreal$PC1,
  covars_park$PC1,
  na.rm = TRUE
)

pc2_lim <- range(
  covars_boreal$PC2,
  covars_park$PC2,
  na.rm = TRUE
)

###Predict interaction plot function
plot_interaction_stan <- function(
    pred_df,
    region,
    ylab = "Abundance",
    title = NULL,
    ylim = NULL,
    xlim = NULL,
    show_ylab = TRUE){
  
  library(colorspace)
  
  region_colors <- list(
    boreal = "#1b9e77",
    park = "#d95f02"
  )
  
  base_col <- region_colors[[region]]
  
  col_low <- lighten(base_col, .4)
  col_high <- darken(base_col, .4)
  
  p <- ggplot(
    pred_df,
    aes(
      x = focal,
      y = pred,
      colour = hf_group,
      fill = hf_group
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      alpha = .10,
      colour = NA
    ) +
    
    geom_line(
      linewidth = 1
    ) +
    
    scale_colour_manual(
      values = c(
        "Low HF" = col_low,
        "High HF" = col_high
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Low HF" = col_low,
        "High HF" = col_high
      )
    ) +
    
    labs(
      x = NULL,
      y = if(show_ylab) ylab else NULL,
      title = title
    ) +
    
    theme_classic() +
    
    theme(
      
      legend.position = "none",
      
      axis.title.y =
        if(show_ylab)
          element_text()
      else
        element_blank(),
      
      axis.text.y =
        if(show_ylab)
          element_text()
      else
        element_blank(),
      
      axis.ticks.y =
        if(show_ylab)
          element_line()
      else
        element_blank()
      
    )
  
  if(!is.null(ylim)){
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  if(!is.null(xlim)){
    p <- p + scale_x_continuous(limits = xlim)
  }
  
  p
}


#Calculate y limits for each species/region to ensure consistent y axis
get_ylim_stan <- function(pc1_preds,
                           pc2_preds){
  
  range(
    
    pc1_preds$lower,
    pc1_preds$upper,
    
    pc2_preds$lower,
    pc2_preds$upper,
    
    na.rm = TRUE
    
  )
  
}

ylim_coyote_boreal <-
  get_ylim_stan(
    coyote_boreal_pc1,
    coyote_boreal_pc2
  )

ylim_coyote_park <-
  get_ylim_stan(
    coyote_park_pc1,
    coyote_park_pc2
  )

ylim_redfox_boreal <-
  get_ylim_stan(
    redfox_boreal_pc1,
    redfox_boreal_pc2
  )

ylim_redfox_park <-
  get_ylim_stan(
    redfox_park_pc1,
    redfox_park_pc2
  )

ylim_lynx_boreal <-
  get_ylim_stan(
    lynx_boreal_pc1,
    lynx_boreal_pc2
  )


#Create Plots for all species/region models
plots <- list(
  
  coyote_boreal_pc1 =
    
    plot_interaction_stan(
      pred_df = coyote_boreal_pc1,
      region = "boreal",
      ylab = "Abundance",
      ylim = ylim_coyote_boreal,
      xlim = pc1_lim,
      show_ylab = TRUE
    ),
  
  coyote_boreal_pc2 =
    
    plot_interaction_stan(
      pred_df = coyote_boreal_pc2,
      region = "boreal",
      ylim = ylim_coyote_boreal,
      xlim = pc2_lim,
      show_ylab = FALSE
    ),
  
  coyote_park_pc1 =
    
    plot_interaction_stan(
      pred_df = coyote_park_pc1,
      region = "park",
      ylab = "Abundance",
      ylim = ylim_coyote_park,
      xlim = pc1_lim,
      show_ylab = TRUE
    ),
  
  coyote_park_pc2 =
    
    plot_interaction_stan(
      pred_df = coyote_park_pc2,
      region = "park",
      ylim = ylim_coyote_park,
      xlim = pc2_lim,
      show_ylab = FALSE
    ),
  
  redfox_boreal_pc1 =
    
    plot_interaction_stan(
      pred_df = redfox_boreal_pc1,
      ylab = "Abundance",
      region = "boreal",
      ylim = ylim_coyote_boreal,
      xlim = pc1_lim,
      show_ylab = TRUE
    ),
  
  redfox_boreal_pc2 =
    
    plot_interaction_stan(
      pred_df = redfox_boreal_pc2,
      region = "boreal",
      ylim = ylim_coyote_boreal,
      xlim = pc2_lim,
      show_ylab = FALSE
    ),
  
  redfox_park_pc1 =
    
    plot_interaction_stan(
      pred_df = redfox_park_pc1,
      region = "park",
      ylab = "Abundance",
      ylim = ylim_redfox_park,
      xlim = pc1_lim,
      show_ylab = TRUE
    ),
 
   redfox_park_pc2 =
    
    plot_interaction_stan(
      pred_df = redfox_park_pc2,
      region = "park",
      ylim = ylim_redfox_park,
      xlim = pc2_lim,
      show_ylab = FALSE
    ),
  
  lynx_boreal_pc1 =
    
    plot_interaction_stan(
      pred_df = lynx_boreal_pc1,
      region = "boreal",
      ylab = "Abundance",
      ylim = ylim_lynx_boreal,
      xlim = pc1_lim,
      show_ylab = TRUE
    ),
  
  lynx_boreal_pc2 =
    
    plot_interaction_stan(
      pred_df = lynx_boreal_pc2,
      region = "boreal",
      ylim = ylim_lynx_boreal,
      xlim = pc2_lim,
      show_ylab = FALSE
    )
)


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




#Create a left side label strip for each row, combine icon and label vertically

make_row_label <- function(species_name,
                           region_name,
                           panel_letter) {
  
  region_col <- case_when(
    region_name == "Boreal/Shield" ~ "#1b9e77",
    TRUE ~ "#d95f02"
  )
  
  icon_path <- file.path(
    "data/phylopics",
    paste0(species_name, ".png")
  )
  
  img <- image_read(icon_path)
  img <- image_colorize(
    img,
    opacity = 100,
    color = region_col
  )
  
  icon_panel <- ggplot() +
    annotation_custom(
      rasterGrob(img),
      xmin = .1,
      xmax = .9,
      ymin = .25,
      ymax = 1
    ) +
    coord_cartesian(
      xlim = c(0,1),
      ylim = c(0,1),
      clip = "off"
    ) +
    theme_void()
  
  text_panel <- ggplot() +
    annotate(
      "text",
      x = 0.05,
      y = 0.95,
      label = panel_letter,
      fontface = "bold",
      size = 6,
      hjust = 0
    ) +
    annotate(
      "text",
      x = .5,
      y = .65,
      label = species_name,
      fontface = "bold",
      size = 4.5
    ) +
    annotate(
      "text",
      x = .5,
      y = .35,
      label = region_name,
      color = region_col,
      size = 3.8
    ) +
    coord_cartesian(
      xlim = c(0,1),
      ylim = c(0,1)
    ) +
    theme_void()
  
  wrap_elements(
    full = icon_panel / text_panel,
    clip = FALSE
  )
}

#Combine row label + PC1 + PC2 into a row
row_coyote_park <- 
  make_row_label("Coyote", "Grassland/Parkland", "A") |
  plots$coyote_park_pc1 |
  plots$coyote_park_pc2

row_coyote_boreal <- 
  make_row_label("Coyote", "Boreal/Shield", "B") |
  plots$coyote_boreal_pc1 |
  plots$coyote_boreal_pc2

row_redfox_park <- 
  make_row_label("Red fox", "Grassland/Parkland", "C") | 
  plots$redfox_park_pc1 | 
  plots$redfox_park_pc2

row_redfox_boreal <- 
  make_row_label("Red fox", "Boreal/Shield", "D") | 
  plots$redfox_boreal_pc1 | 
  plots$redfox_boreal_pc2

row_lynx_boreal <- 
  make_row_label("Canada lynx",   "Boreal/Shield", "E") | 
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

#Build region hf legend
legend_df <- tibble::tribble(
  ~region,     ~x_region, ~x1, ~x2,
  "Boreal",      0.5,      2.0, 3.8,
  "Parkland",    0.5,      2.0, 3.8
)

make_legend_row <- function(region_name,
                            low_col,
                            high_col) {
  
  ggplot() +
    
    # Region label
    annotate(
      "text",
      x = 0.3,
      y = 1,
      label = region_name,
      hjust = 0,
      fontface = "bold",
      size = 3.8
    ) +
    
    # Low HF
    annotate(
      "point",
      x = 2.0,
      y = 1,
      colour = low_col,
      size = 3.5
    ) +
    
    annotate(
      "text",
      x = 2.25,
      y = 1,
      label = "Low HF",
      hjust = 0,
      size = 3.2
    ) +
    
    # High HF
    annotate(
      "point",
      x = 3.8,
      y = 1,
      colour = high_col,
      size = 3.5
    ) +
    
    annotate(
      "text",
      x = 4.05,
      y = 1,
      label = "High HF",
      hjust = 0,
      size = 3.2
    ) +
    
    coord_cartesian(
      xlim = c(0, 5.5),
      ylim = c(0.8, 1.2),
      clip = "off"
    ) +
    
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
}


legend_boreal <- make_legend_row(
  "Boreal",
  lighten("#1b9e77", .4),
  darken("#1b9e77", .4)
)

legend_park <- make_legend_row(
  "Grassland",
  lighten("#d95f02", .4),
  darken("#d95f02", .4)
)

hf_region_legend <-
  legend_boreal /
  legend_park +
  plot_layout(
    heights = c(1,1)
  )
#Combine grid + legend with controlled widths

top_row <-
  hf_region_legend |
  plot_spacer()


final <-
  top_row /
  grid /
  (plot_spacer() | xlab_pc1 | xlab_pc2) +
  plot_layout(
    widths = c(0.6, 3, 3),
    heights = c(0.05, 1, 0.04)
  )

final

ggsave("interactionplots_abmi.pdf", final, width = 7, height = 12)
ggsave("interactionplots_abmi.png", final, width = 7, height = 12)

