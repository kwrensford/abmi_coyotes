#Building forest plots
library(vroom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ubms)
library(tidyverse)
library(ggimage)
library(patchwork)
library(cowplot)
library(magick)

coyote_nmix_ubms_boreal <- readRDS("coyote_nmix_ubms_boreal.rds")
coyote_nmix_ubms_park <- readRDS("coyote_nmix_ubms_park.rds")
redfox_nmix_stan_boreal <- readRDS("redfox_nmix_stan_boreal.rds")
redfox_nmix_stan_park <- readRDS("redfox_nmix_stan_park.rds")
lynx_nmix_stan <- readRDS("lynx_nmix_stan.rds")

coyote_boreal_est <- summary(coyote_nmix_ubms_boreal, submodel = "state")
coyote_park_est <- summary(coyote_nmix_ubms_park, submodel = "state")
redfox_boreal_est <- summary(redfox_nmix_stan_boreal, submodel = "state")
redfox_park_est <- summary(redfox_nmix_stan_park, submodel = "state")
lynx_boreal_est <- summary(lynx_nmix_stan, submodel = "state")

extract_estimates <- function(df, species, region){
  
  df %>%
    rownames_to_column("Variable") %>%
    filter(Variable != "(Intercept)") %>%   # remove intercept
    transmute(
      species,
      region,
      Variable,
      mean = mean,
      lower = `2.5%`,
      upper = `97.5%`
    )
  
}

all_estimates <- bind_rows(
  
  extract_estimates(coyote_boreal_est, "Coyote", "Boreal"),
  extract_estimates(coyote_park_est,   "Coyote", "Parkland"),
  
  extract_estimates(redfox_boreal_est, "Red Fox", "Boreal"),
  extract_estimates(redfox_park_est,   "Red Fox", "Parkland"),
  
  extract_estimates(lynx_boreal_est,   "Lynx", "Boreal")
  
)

all_estimates <- all_estimates %>%
  mutate(
    excludes_zero = !(lower < 0 & upper > 0)
  )

all_estimates <- all_estimates %>%
  mutate(
    Variable = factor(
      Variable,
      levels = rev(c(
        "PC1",
        "PC2",
        "hf",
        "PC1:hf",
        "hf:PC2"
      ))
    )
  )

region_colors <- list(
  boreal = "#1b9e77",
  park   = "#d95f02"
)

#Species images
icon_df <- tibble(
  species = c(
    "Coyote",
    "Red Fox",
    "Lynx"
  ),
  image = c(
    "data/phylopics/Coyote.png",
    "data/phylopics/Red Fox.png",
    "data/phylopics/Canada lynx.png"
  ),
  x = 1:3,
  y = 1
)

icon_plot <- ggplot(
  icon_df,
  aes(x, y)
) +
  geom_image(
    aes(image = image),
    size = 0.12
  ) +
  geom_text(
    aes(label = species),
    y = 0.72,
    fontface = "bold",
    size = 4
  ) +
  theme_void() +
  coord_cartesian(
    xlim = c(0.5, 3.5),
    ylim = c(0.5, 1.2)
  )

forest_plot <- ggplot(
  all_estimates,
  aes(
    x = mean,
    y = Variable,
    colour = region
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey50"
  ) +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0,
    position = position_dodge(width = 0.6)
  ) +
  
  geom_point(
    size = 3,
    position = position_dodge(width = 0.6)
  ) +
  
  facet_wrap(
    ~ species,
    nrow = 1
  ) +
  
  scale_colour_manual(
    values = c(
      Boreal = "#1b9e77",
      Parkland = "#d95f02"
    )
  ) +
  
  theme_classic() +
  
  labs(
    x = "Posterior estimate (95% CI)",
    y = NULL
  ) 

forest_plot

ggdraw(forest_plot) +
  draw_image("data/phylopics/Coyote.png",
             x = 0.12, y = 0.92,
             width = 0.08) +
  draw_image("data/phylopics/Red fox.png",
             x = 0.45, y = 0.92,
             width = 0.08) +
  draw_image("data/phylopics/Canada lynx.png",
             x = 0.78, y = 0.92,
             width = 0.08)

final_plot <-
  icon_plot /
  forest_plot +
  plot_layout(
    heights = c(0.12, 1)
  )

final_plot
