#CLimate Variable Dimension Reduction with Principal Component Analysis

#Load in covariate data
##Filter out climate covariates

covariates_climate <- covariates_unique %>%
  dplyr::select(location, Latitude, Longitude, AHM, PET, FFP, MAP, MAT, MWMT, MCMT, region)

##Keep location info
covariates_climate <- covariates_climate %>%
  filter(region != "Rocky+Foothills")

covar_locations <- covariates_climate %>%
  dplyr::select(location, Latitude, Longitude, region)

##Filter to climate variables only
climate_var <- covariates_climate %>%
  dplyr::select(-location, -Latitude, -Longitude, -region)

#Scale climate variabes
climate_scaled <- scale(climate_var)

##Run PCA
pca_res <- prcomp(climate_scaled, center = TRUE, scale. = TRUE)

plot(pca_res, type = "lines")

summary(pca_res)

biplot(pca_res)


##Extract principal components
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
  select(x, y) %>% 
  bind_cols(climate_sites)

# Extract scaled loadings 
# autoplot always puts loadings in the LAST layer
loadings_df <- pb$data[[length(pb$data)]]

# Clean up loadings_df
loadings_df <- loadings_df %>%
  select(x, y, label)

# Clean up loadings_df
loadings_df <- loadings_df %>%
  select(x, y, label)

#Compute complex hulls for region polygons
hulls <- climate_sites %>%
  group_by(region) %>%
  slice(chull(PC1, PC2))   # convex hull for each region

hulls <- plot_data %>%
  group_by(region) %>%
  slice(chull(x, y))

#Build final plot with region hulls and arrow labels
pc_plot <- pc_plot + 
  geom_polygon(
    data = hulls,
    aes(x = x, y = y, group = region, linetype = region),
    fill = NA,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = loadings_df,
    aes(x = x, y = y, label = label),
    fontface = "bold",
    size = 4,
    color = "black",
    box.padding = 0,
    point.padding = 0,
    min.segment.length = 0,
    force = 0.0001,
    max.time = 0.1,
    segment.color = "grey40",
    max.overlaps = Inf
  )

pc_plot

ggsave("data/pc_plot.png", pc_plot, width = 8, height = 7)