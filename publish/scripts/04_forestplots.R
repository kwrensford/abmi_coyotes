#Generate Forest Plots

library(tidyverse)

extract_posteriors <- function(file,
                               species,
                               region){
  
  read_csv(file, show_col_types = FALSE) %>%
    
    filter(str_detect(variable, "^beta_")) %>%
    
    transmute(
      species,
      region,
      
      Variable = recode(
        variable,
        
        beta_clim1     = "PC1",
        beta_clim2     = "PC2",
        beta_human     = "hf",
        beta_interact1 = "PC1:hf",
        beta_interact2 = "hf:PC2"
      ),
      
      mean  = mean,
      lower = q5,
      upper = q95
    )
}

all_estimates <- bind_rows(
  
  extract_posteriors(
    "output/summaries/coyote_boreal_summary.csv",
    "Coyote",
    "Boreal"
  ),
  
  extract_posteriors(
    "output/summaries/coyote_park_summary.csv",
    "Coyote",
    "Parkland"
  ),
  
  extract_posteriors(
    "output/summaries/redfox_boreal_summary.csv",
    "Red Fox",
    "Boreal"
  ),
  
  extract_posteriors(
    "output/summaries/redfox_park_summary.csv",
    "Red Fox",
    "Parkland"
  ),
  
  extract_posteriors(
    "output/summaries/lynx_boreal_summary.csv",
    "Lynx",
    "Boreal"
  )
)

#Save as csv

write.csv(all_estimates, file = "data/all_estimates_stan.csv")

#Flag significant effects
all_estimates <- all_estimates %>%
  
  mutate(
    excludes_zero =
      !(lower < 0 & upper > 0)
  )

#Order variables
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

#Remove intercepts
all_estimates <- all_estimates %>%
  filter(Variable != "NA")

#Set dodge position for points
pd <- position_dodge(width = 0.6)



#Forest plot faceted by region

forest_plot <- ggplot(
  all_estimates,
  aes(
    x = mean,
    y = Variable,
    colour = species,
    group = species
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey50"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    position = pd,
    height = 0
  ) +
  
  geom_point(
    aes(shape = species),
    position = pd,
    size = 3
  ) +
  
  facet_wrap(
    ~ region,
    nrow = 1
  ) +
  
  labs(
    x = "Posterior estimate",
    y = NULL
  ) +
  
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.y.left = element_text(size = 16),
    strip.text = element_text(size = 16))

forest_plot

ggsave("forest_plot.png", forest_plot)
ggsave("forest_plot.pdf", forest_plot)
