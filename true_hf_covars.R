#Download and organiz True HF Data

library(dplyr)
library(tidyr)
library(ggplot2)
library(vroom)

hf_true <- vroom("data/chfi_truelocals.csv")

hf_public <- vroom("data/chfi_publiclocals.csv")

covars <- vroom("data/covars.csv")

plot(covars$hf, hf_true$hf)
plot(covars$hf, hf_public$hf)
plot(hf_public$hf, hf_true$hf)

summary(lm(hf_true$hf~covars$hf))

#Ensure orders match
hf_public <- hf_public[match(covars$location, hf_public$location), ]
hf_true <- hf_true[match(covars$location, hf_true$location), ]

#Create covariate object for public hf data
covars_public <- covars

#Create covariate object for true hf data
covars <- covars %>%
  dplyr::select(location, site, cam_orientation, Latitude, Longitude, region, PC1, PC2) %>%
  left_join(dplyr::select(hf_true, location, hf), by = "location")

covars <- covars %>%
  mutate(
    baited = cam_orientation %in% c("NE", "SW")
  )

covars_boreal <- covars %>%
  filter(region == "Shield/Boreal")

covars_park <- covars %>%
  filter(region == "Grassland/Parkland")

closeAllConnections()

gc()

write.csv(covars, file = "data/covars.csv", row.names = FALSE)
write.csv(covars_public, file = "data/covars_public.csv", row.names = FALSE)
write.csv(covars_boreal, file = "data/covars_boreal.csv", row.names = FALSE)
write.csv(covars_park, file = "data/covars_park.csv", row.names = FALSE)

