##Human footprint covariate extraction for Marcus

library(vroom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(forcats)
library(lubridate)
library(unmarked)
library(camtrapR)
library(overlap)
library(ggmap)
library(reshape2)
library(corrplot)
library(terra)
library(mapview)
library(sf)
library(exactextractr)
library(xml2)
library(geodata)
library(viridis)
library(scales)
library(colorspace)

#Load in site covar data (adapt if you use a different file, rename)
covars <- read.csv("data/covars.csv")
covar_locations <- covars %>%
  select(location, site, Latitude, Longitude)

# Load rasters (human footprint TIFF) and crop to alberta provincial boundary
hf_raster <- rast("data/chfi/cum_threat2020.02.18.tif")
roads_raster <- rast("data/chfi/roads.tif")
rail_raster <- rast("data/chfi/rail.tif")
built_raster <- rast("data/chfi/built.tif")
oil_gas_raster <- rast("data/chfi/oil_gas.tif")
forestry_raster <- rast("data/chfi/forestry_harvest.tif")
mines_raster <- rast("data/chfi/mines.tif")
crop_raster <- rast("data/chfi/crop.tif")
pasture_raster <- rast("data/chfi/pasture.tif")

roads_raster <- resample(roads_raster, hf_raster)
rail_raster <- resample(rail_raster, hf_raster)
built_raster <- resample(built_raster, hf_raster)
oil_gas_raster <- resample(oil_gas_raster, hf_raster)
forestry_raster <- resample(forestry_raster, hf_raster)
mines_raster <- resample(mines_raster, hf_raster)
crop_raster <- resample(crop_raster, hf_raster)
pasture_raster <- resample(pasture_raster, hf_raster)

hf_stack <- c(hf_raster, roads_raster, rail_raster, built_raster, oil_gas_raster, forestry_raster, mines_raster, crop_raster, pasture_raster)

names(hf_stack) <- c("hf", "roads", "rail", "built", "oil_gas", "forestry", "mines", "crop", "pasture")

##Obtain Canada and Alberta boundary shapefiles

url <- "https://geospatial.alberta.ca/titan/rest/services/boundary/goa_administrative_area/MapServer/0"

canada <- gadm("Canada", level = 1, path = "data/")
alberta <- canada[canada$NAME_1 == "Alberta", ]
alberta_proj <- project(alberta, hf_stack)

hf_alberta <- crop(hf_stack, alberta_proj)
hf_alberta_masked <- mask(hf_alberta, alberta_proj)

# Convert sites to spatial object
covar_locations <- covar_locations %>%
  filter(!is.na(Longitude) & !is.na(Latitude))

sites_hf <- st_as_sf(covar_locations, coords = c("Longitude", "Latitude"), crs = 4326)
sites_hf <- st_transform(sites_hf, crs(hf_alberta_masked))
sites_vect <- vect(sites_hf)

plot(hf_alberta_masked$hf)
plot(sites_vect, add = TRUE, col = "darkorange", cex = 0.5)

#Define buffer size
buffers <- st_buffer(sites_hf, dist = 5000)
buffers_vect <- vect(buffers)

#Linear Features
hf_linear <- hf_alberta_masked$rail + hf_alberta_masked$roads
names(hf_linear) <- "linear"

plot(hf_linear)
plot(sites_vect, add = TRUE, col = "darkorange", cex = 0.5)

#Industrial
hf_industrial <- hf_alberta_masked$oil_gas + hf_alberta_masked$mines + hf_alberta_masked$forestry
names(hf_industrial) <- "industrial"

plot(hf_industrial)
plot(sites_vect, add = TRUE, col = "darkorange", cex = 0.5)

#Agriculture
hf_agriculture <- hf_alberta_masked$crop + hf_alberta_masked$pasture
names(hf_agriculture) <- "agriculture"

plot(hf_agriculture)
plot(sites_vect, add = TRUE, col = "darkorange", cex = 0.5)

#Extract data
hf_vals <- terra::extract(hf_alberta_masked, buffers_vect, fun = mean, na.rm = TRUE)
linear_vals <- terra::extract(hf_linear, buffers_vect, fun = mean, na.rm = TRUE)
industrial_vals <- terra::extract(hf_industrial, buffers_vect, fun = mean, na.rm = TRUE)
agriculture_vals <- terra::extract(hf_agriculture, buffers_vect, fun = mean, na.rm = TRUE)

#Combine hf data with 
covar_locations_hf <- bind_cols(covar_locations, hf_vals[,-1])

write.csv(covar_locations_hf, file = "chfi_truelocals.csv", row.names = FALSE)
