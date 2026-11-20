#PCA Plots

library(tidyverse)
library(ggfortify)
library(gggrepel)
library(ggnewscale)
        
covariates_climate <- read.csv("data/pcacovars.csv")

covariates_climate <- covariates_climate %>%
  mutate(region = case_when(
    region == "Shield/Boreal" ~ "Boreal/Shield",
    TRUE ~ region
  ))
        
covar_locations <- covariates_climate %>%
  dplyr::select(location, Latitude, Longitude, region)

##Climate variables only
climate_var <- covariates_climate %>%
    dplyr::select(-location, -Latitude, -Longitude, -region)

climate_scaled <- scale(climate_var)

##Run PCA
pca_res <- prcomp(climate_scaled, center = TRUE, scale. = TRUE)

plot(pca_res, type = "lines")

summary(pca_res)

biplot(pca_res)


##Extract PC's
climate_pcs <- as.data.frame(pca_res$x[, 1:2])
colnames(climate_pcs) <- c("PC1", "PC2")

# Loadings (variable contributions to PCs)
round(pca_res$rotation[, 1:2], 3)
load_df <- as.data.frame(pca_res$rotation[, 1:2])
load_df$label <- rownames(load_df)

# Bind with site info
climate_sites <- cbind(covar_locations, climate_pcs)

#Plot PC's
mult <- min(
  (max(autoplot(pca_res)$data$PC1) - min(autoplot(pca_res)$data$PC1)) /
    (max(load_df$PC1) - min(load_df$PC1)),
  (max(autoplot(pca_res)$data$PC2) - min(autoplot(pca_res)$data$PC2)) /
    (max(load_df$PC2) - min(load_df$PC2))
)

load_df[, c("PC1", "PC2")] <- load_df[, c("PC1", "PC2")] * mult


pc_plot <- autoplot(pca_res, data = climate_sites, color = "Latitude",
         loadings = TRUE, loadings.color = "black", 
         loadings.label = FALSE)+
  geom_point(aes(color = Latitude), size = 3) +
  scale_color_viridis_c(option = "plasma", direction = -1)+
  labs(x = "PC1 - Overall Climate",
       y = "PC2 - Aridity")+
  theme_classic() +
  theme(axis.title = element_text(size = 15))

plot_data <- ggplot_build(pc_plot)$data[[1]]

plot_data <- bind_cols(
  climate_sites, 
  plot_data[,c("x","y")] #autoplot-scaled PC1/PC2
)

#Extract scaled PCA scores
pb <- ggplot_build(pc_plot)

scores_df <- pb$data[[1]] %>% 
  dplyr::select(x, y) %>% 
  bind_cols(climate_sites)


loadings_df <- load_df %>%
  dplyr::rename(
    x = PC1,
    y = PC2
  )

#Compute complex hulls for region polygons
hulls <- climate_sites %>%
  group_by(region) %>%
  slice(chull(PC1, PC2))   # convex hull for each region

hulls <- plot_data %>%
  group_by(region) %>%
  slice(chull(x, y))

hulls <- hulls %>%
  mutate(
    hull_col = case_when(
      region == "Boreal/Shield" ~ "#1b9e77",
      region == "Grassland/Parkland" ~ "#d95f02"
    )
  )

hulls$region <- factor(
  hulls$region,
  levels = c("Boreal/Shield", "Grassland/Parkland")
)


#Build final plot with region hulls and arrow labels
pc_plot <- pc_plot + 
  
  new_scale_color() +
  
  geom_polygon(
    data = hulls,
    aes(
      x = x,
      y = y,
      group = region,
      color = region,
      linetype = region
    ),
    fill = NA,
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    name = "Region",
    values = c(
      "Boreal/Shield" = "#1b9e77",
      "Grassland/Parkland" = "#d95f02"
    )
  ) +scale_linetype_manual(
    name = "Region",
    values = c(
      "Boreal/Shield" = "solid",
      "Grassland/Parkland" = "longdash"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c("solid", "longdash"),
        linewidth = 1
      )
    ),
    linetype = "none"
  ) +
  geom_text(
    data = loadings_df,
    aes(
      x = x ,
      y = y ,
      label = label
    ),
    fontface = "bold",
    size = 4
  )

pc_plot

ggsave("data/pc_plot.png", pc_plot, width = 8, height = 7)
