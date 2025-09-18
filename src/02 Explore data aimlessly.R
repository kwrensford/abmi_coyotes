library(vroom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(lubridate)
library(unmarked)
library(camtrapR)
library(overlap)
library(ggmap)
library(reshape2)

# Camera Reports
all_main_reports <- vroom::vroom("data/ABMI EH 2014-2019 Main Reports.csv")
all_image_reports <- vroom::vroom("data/ABMI EH 2014-2019 Image Reports.csv")

head(all_main_reports)
head(all_image_reports)

camlocations <- all_main_reports %>%
  distinct(location, .keep_all = TRUE)

##Create column for year of report
all_main_reports$year <- lubridate::year(all_main_reports$image_date_time)
all_image_reports$year <- lubridate::year(all_image_reports$image_date_time)

##Filter reports to 2015-2019 range
study_image_reports<-all_image_reports %>%
  filter(year > 2014 & year < 2020)

study_main_reports<-all_main_reports %>%
  filter(year > 2014 & year < 2020)

sort(unique(all_main_reports$species_common_name))

##Segment Data for each Year
main_reports_2015<-all_main_reports%>%
  filter(year==2015)

main_reports_2016<-all_main_reports%>%
  filter(year==2016)

main_reports_2017<-all_main_reports%>%
  filter(year==2017)

main_reports_2018<-all_main_reports%>%
  filter(year==2018)

main_reports_2019<-all_main_reports%>%
  filter(year==2019)

#Camera effort data
total_operation_days<- vroom::vroom("data/Total Days of Operations ABMI EH 2014-2020.csv")

camera_operation_dates <- vroom::vroom("data/Operation Dates ABmi EH 2014-2018.csv")

##Use image reports file to extract camera operation dates for each location
camera_operation_dates_v2<-study_image_reports %>%
  group_by(year, project_id, location) %>%
  summarize(
    start_date_time = min(image_date_time),
    end_date_time = max(image_date_time)
  )

camera_operation_dates_v2$start_date_time<-lubridate::as_date(camera_operation_dates_v2$start_date_time)
camera_operation_dates_v2$end_date_time<-lubridate::as_date(camera_operation_dates_v2$end_date_time)

##What threshold of total operation dates should we exclude cameras from analyses?

##Does the dataframe for an unmarked analysis need to be symmetrical?

#Other attempt to make detection matrix for unmarked

# Naive occupancy

##  Coyotes
coyote_site_occ <- study_main_reports %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote <- mean(coyote_site_occ$coyotedetection)

ggplot(coyote_site_occ, aes(x = longitude, y = latitude, color = factor(coyotedetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Coyote Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Coyote Detection")

##  Gray Wolf
wolf_site_occ <- study_main_reports %>%
  group_by(location, longitude, latitude)%>%
  summarise(wolfdetection = as.integer(any(species_common_name == "Gray Wolf")))

naive_occ_wolf <- mean(wolf_site_occ$wolfdetection)

ggplot(wolf_site_occ, aes(x = longitude, y = latitude, color = factor(wolfdetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Wolf Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Wolf Detection")

##  Lynx
lynx_site_occ <- study_main_reports %>%
  group_by(location, longitude, latitude)%>%
  summarise(lynxdetection = as.integer(any(species_common_name == "Canada Lynx")))

naive_occ_lynx <- mean(lynx_site_occ$lynxdetection)

ggplot(lynx_site_occ, aes(x = longitude, y = latitude, color = factor(lynxdetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Lynx Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Lynx Detection")

##  Marten
marten_site_occ <- study_main_reports %>%
  group_by(location, longitude, latitude)%>%
  summarise(martendetection = as.integer(any(species_common_name == "Marten")))

naive_occ_marten <- mean(marten_site_occ$martendetection)

ggplot(marten_site_occ, aes(x = longitude, y = latitude, color = factor(martendetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Marten Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Marten Detection")

##  Fisher
fisher_site_occ <- study_main_reports %>%
  group_by(location, longitude, latitude)%>%
  summarise(fisherdetection = as.integer(any(species_common_name == "Fisher")))

naive_occ_fisher <- mean(fisher_site_occ$fisherdetection)

ggplot(fisher_site_occ, aes(x = longitude, y = latitude, color = factor(fisherdetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Fisher Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Fisher Detection")

##Combine into a single dataset
naive_detections<-coyote_site_occ %>%
  inner_join(wolf_site_occ)%>%
  inner_join(lynx_site_occ)%>%
  inner_join(marten_site_occ)%>%
  inner_join(fisher_site_occ)

naive_occ <- naive_detections%>%
  ungroup()%>%
  dplyr::select(coyotedetection:fisherdetection)%>%
  colMeans()

print(naive_occ)

##Compute co-occurence

##Coyotelynx overlap
coyotelynx_overlap<-naive_detections%>%
  dplyr::select(location, longitude, latitude, coyotedetection, lynxdetection)

coyotelynx_overlap <- coyotelynx_overlap %>%
  mutate(OverlapCategory = case_when(
    coyotedetection == 1 & lynxdetection == 1 ~ "Both Detected",
    coyotedetection == 1 & lynxdetection == 0 ~ "Coyote Only",
    coyotedetection == 0 & lynxdetection == 1 ~ "Lynx Only",
    TRUE ~ "Neither"  # Optional: remove this line if you only want sites with detections
  ))

# Percentage of site where species cooccur
coyotelynx_naiveocc <- sum(coyotelynx_overlap$OverlapCategory=="Both Detected")/length(coyotelynx_overlap$location)

overlap_colors_lynx <- c("Both Detected" = "purple", "Coyote Only" = "red", "Lynx Only" = "blue", "Neither" = "gray")

ggplot(coyotelynx_overlap, aes(x = longitude, y = latitude, color = OverlapCategory)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = overlap_colors_lynx) +
  theme_minimal() +
  labs(title = "Coyote-Lynx Co-Occurrence Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Overlap Category")

##Coyotewolf overlap
coyotewolf_overlap<-naive_detections%>%
  dplyr::select(location, longitude, latitude, coyotedetection, wolfdetection)

coyotewolf_overlap <- coyotewolf_overlap %>%
  mutate(OverlapCategory = case_when(
    coyotedetection == 1 & wolfdetection == 1 ~ "Both Detected",
    coyotedetection == 1 & wolfdetection == 0 ~ "Coyote Only",
    coyotedetection == 0 & wolfdetection == 1 ~ "Wolf Only",
    TRUE ~ "Neither"  # Optional: remove this line if you only want sites with detections
  ))

coyotewolf_naiveocc <- sum(coyotewolf_overlap$OverlapCategory=="Both Detected")/length(coyotewolf_overlap$location)

overlap_colors_wolf <- c("Both Detected" = "purple", "Coyote Only" = "red", "Wolf Only" = "blue", "Neither" = "gray")

ggplot(coyotewolf_overlap, aes(x = longitude, y = latitude, color = OverlapCategory)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = overlap_colors_wolf) +
  theme_minimal() +
  labs(title = "Coyote-Wolf Co-Occurrence Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Overlap Category")

##Coyotemarten
coyotemarten_overlap<-naive_detections%>%
  dplyr::select(location, longitude, latitude, coyotedetection, martendetection)

coyotemarten_overlap <- coyotemarten_overlap %>%
  mutate(OverlapCategory = case_when(
    coyotedetection == 1 & martendetection == 1 ~ "Both Detected",
    coyotedetection == 1 & martendetection == 0 ~ "Coyote Only",
    coyotedetection == 0 & martendetection == 1 ~ "Marten Only",
    TRUE ~ "Neither"  # Optional: remove this line if you only want sites with detections
  ))

coyotemarten_naiveocc <- sum(coyotemarten_overlap$OverlapCategory=="Both Detected")/length(coyotemarten_overlap$location)

overlap_colors_marten <- c("Both Detected" = "purple", "Coyote Only" = "red", "Marten Only" = "blue", "Neither" = "gray")

ggplot(coyotemarten_overlap, aes(x = longitude, y = latitude, color = OverlapCategory)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = overlap_colors_marten) +
  theme_minimal() +
  labs(title = "Coyote-Marten Co-Occurrence Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Overlap Category")

##Coyotefisher
coyotefisher_overlap<-naive_detections%>%
  dplyr::select(location, longitude, latitude, coyotedetection, fisherdetection)

coyotefisher_overlap <- coyotefisher_overlap %>%
  mutate(OverlapCategory = case_when(
    coyotedetection == 1 & fisherdetection == 1 ~ "Both Detected",
    coyotedetection == 1 & fisherdetection == 0 ~ "Coyote Only",
    coyotedetection == 0 & fisherdetection == 1 ~ "Fisher Only",
    TRUE ~ "Neither"  # Optional: remove this line if you only want sites with detections
  ))

coyotefisher_naiveocc <- sum(coyotefisher_overlap$OverlapCategory=="Both Detected")/length(coyotefisher_overlap$location)

overlap_colors_fisher <- c("Both Detected" = "purple", "Coyote Only" = "red", "Fisher Only" = "blue", "Neither" = "gray")

ggplot(coyotefisher_overlap, aes(x = longitude, y = latitude, color = OverlapCategory)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = overlap_colors_fisher) +
  theme_minimal() +
  labs(title = "Coyote-Fisher Co-Occurrence Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Overlap Category")

## Partition data by year and plot.
### Coyote
coyote_site_occ_2015 <- main_reports_2015 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2015 <- mean(coyote_site_occ_2015$coyotedetection)

ggplot(coyote_site_occ_2015, aes(x = longitude, y = latitude, color = factor(coyotedetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Coyote Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Coyote Detection")

coyote_site_occ_2016 <- main_reports_2016 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2016 <- mean(coyote_site_occ_2016$coyotedetection)

ggplot(coyote_site_occ_2016, aes(x = longitude, y = latitude, color = factor(coyotedetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Coyote Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Coyote Detection")


coyote_site_occ_2017 <- main_reports_2017 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2017 <- mean(coyote_site_occ_2017$coyotedetection)

ggplot(coyote_site_occ_2017, aes(x = longitude, y = latitude, color = factor(coyotedetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Coyote Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Coyote Detection")


coyote_site_occ_2018 <- main_reports_2018 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2018 <- mean(coyote_site_occ_2018$coyotedetection)

ggplot(coyote_site_occ_2018, aes(x = longitude, y = latitude, color = factor(coyotedetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Coyote Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Coyote Detection")


coyote_site_occ_2019 <- main_reports_2019 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2019 <- mean(coyote_site_occ_2019$coyotedetection)

ggplot(coyote_site_occ_2019, aes(x = longitude, y = latitude, color = factor(coyotedetection))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray", "1" = "red"), labels = c("Not Detected", "Detected")) +
  theme_minimal() +
  labs(title = "Coyote Detection Across Sites",
       x = "Longitude", y = "Latitude",
       color = "Coyote Detection")


### Wolf
wolf_site_occ_2015 <- main_reports_2015 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_wolf_2015 <- mean(wolf_site_occ_2015$coyotedetection)

wolf_site_occ_2016 <- main_reports_2016 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2016 <- mean(coyote_site_occ_2016$coyotedetection)

coyote_site_occ_2017 <- main_reports_2017 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2017 <- mean(coyote_site_occ_2017$coyotedetection)

coyote_site_occ_2018 <- main_reports_2018 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2018 <- mean(coyote_site_occ_2018$coyotedetection)

coyote_site_occ_2019 <- main_reports_2019 %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote_2019 <- mean(coyote_site_occ_2019$coyotedetection)
### Lynx
### Marten
### Fisher

#Covariate exploration

## Human Footprint
ehs_humanfootprint<-vroom::vroom("data/Point Level Vegetation HF for Kwasi.csv")
ehs_humanfootprint_buffered<-vroom::vroom("data/150m Buffer Vegetation HF for Kwasi.csv")

buffered_renamed <- ehs_humanfootprint_buffered %>%
  rename_with(~ paste0(.x, "_buf"),
              .cols = -location)

ehs_humanfootprint_merged <- ehs_humanfootprint %>%
  left_join(buffered_renamed, by = "location")

## Climate variables
#1) Annual variables:

#Directly calculated annual variables:

#MAT              mean annual temperature (°C),

#MWMT           mean warmest month temperature (°C),

#MCMT            mean coldest month temperature (°C),

#TD                   temperature difference between MWMT and MCMT, or continentality (°C),

#MAP               mean annual precipitation (mm),

#MSP                May to September precipitation (mm),

#AHM  annual heat-moisture index (MAT+10)/(MAP/1000))

#SHM               summer heat-moisture index ((MWMT)/(MSP/1000))

#FFP                Frost free period

#PET                Potential evapotranspiration

## Climate Plotting

## MCMT
ggplot(covariates_fixed, aes(x = Longitude, y = Latitude, color = MCMT)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()

ggplot(ehs_humanfootprint, aes(x = Lat, y = MCMT))+
  geom_point()+
  theme_minimal()

## MWMT
ggplot(ehs_humanfootprint, aes(x = Long, y = Lat, color = MWMT)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()

ggplot(ehs_humanfootprint, aes(x = Lat, y = MWMT))+
  geom_point()+
  theme_minimal()

##FFP
ggplot(ehs_humanfootprint, aes(x = Long, y = Lat, color = FFP)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()

ggplot(ehs_humanfootprint, aes(x = Lat, y = FFP))+
  geom_point()+
  theme_minimal()

## Human Footprint Veg plotting
ggplot(ehs_humanfootprint, aes(x = Long, y = Lat, color = FFP)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()

ggplot(ehs_humanfootprint, aes(x = Lat, y = FFP))+
  geom_point()+
  theme_minimal()

## Human footprint
hist(ehs_humanfootprint_merged)

ggplot(ehs_humanfootprint_buffered, aes(x = Long, y = Lat, color = RoadHardSurface)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()

ggplot(ehs_humanfootprint_merged, aes(x = Lat, y = RoadHardSurface_buf))+
  geom_point()+
  theme_minimal()

## Human footprint x climate
ggplot(ehs_humanfootprint_merged, aes(x = MWMT, y = RoadHardSurface_buf ))+
  geom_point()+
  theme_minimal()

ggplot(ehs_humanfootprint_merged, aes(x = MCMT, y = RoadHardSurface_buf ))+
  geom_point()+
  theme_minimal()

ggplot(ehs_humanfootprint_merged, aes(x = MCMT, y = RuralResidentialIndustrial ))+
  geom_point()+
  theme_minimal()

## Covariate relationships with naive occupancy

## ClimateXcoyote
naive_data <- naive_detections %>%
  left_join(covariates_fixed, by = "location")

ggplot(naive_data, aes(x = factor(coyotedetection), y = MAT)) +
  geom_boxplot() +
  labs(x = "Coyote Detected (0 = No, 1 = Yes)", y = "Mean Annual Temperature")

ggplot(naive_data, aes(x = factor(coyotedetection), y = MAP)) +
  geom_boxplot() +
  labs(x = "Coyote Detected (0 = No, 1 = Yes)", y = "Mean Annual Precipitation")

ggplot(naive_data, aes(x = factor(coyotedetection), y = FFP)) +
  geom_boxplot() +
  labs(x = "Coyote Detected (0 = No, 1 = Yes)", y = "Frost-Free Period")

ggplot(naive_data, aes(x = factor(coyotedetection), y = MCMT)) +
  geom_boxplot() +
  labs(x = "Coyote Detected (0 = No, 1 = Yes)", y = "Mean Temp Coldest Month")

ggplot(naive_data, aes(x = factor(coyotedetection), y = MWMT)) +
  geom_boxplot() +
  labs(x = "Coyote Detected (0 = No, 1 = Yes)", y = "Mean Temp Warmest Month")

ggplot(naive_data, aes(x = factor(coyotedetection), y = RoadHardSurface_buf)) +
  geom_boxplot() +
  labs(x = "Coyote Detected (0 = No, 1 = Yes)", y = "RoadHardSurface (buffered)")

## Prelim models
naive_model <- glm(coyotedetection ~ MAT*latitude+RoadHardSurface_buf+FFP*latitude,
                   family = binomial, data = naive_data)

summary(naive_model)

##Interaction between MAT and latitude stronger than effect of latitude alone, but effect size small.
##Need to account for elevation, probably explains a lot of FFP's effect

naive_model_lynx <- glm(lynxdetection ~ coyotedetection*FFP + MAT*latitude+FFP*latitude,
                        family = binomial, data = naive_data)

summary(naive_model_lynx)
#Slight interaction of ffp and coyote detection effect on lynx.
#I THINK the effect is that as frost free period increases, effect of coyote detection becomes negative? So coyotes less likely to

naive_model_marten <- glm(martendetection ~ coyotedetection + fisherdetection + MAT*latitude+FFP,
                        family = binomial, data = naive_data)

summary(naive_model_marten)

##Marten only species where coyote detection significantly predicts marten detection, v interesting

naive_model_fisher <- glm(fisherdetection ~ coyotedetection + MAT*latitude,
                          family = binomial, data = naive_data)

summary(naive_model_fisher)

##Buffered dataset
naive_model_buffered <- glm(coyotedetection ~ latitude+UrbanIndustrial+UrbanResidence+RoadHardSurface,
                            family = binomial, data = naive_data_buffered)

summary(naive_model_buffered)

## RoadHardSurface predicts coyote occupancy, but only with 150m buffered data

naive_model_marten_buffered <- glm(martendetection ~ latitude+fisherdetection+coyotedetection*RoadHardSurface,
                                   family = binomial, data = naive_data_buffered)

summary(naive_model_marten_buffered)

##Roads important for martens too, but relationshop is weak, coyote presence still strongest predictor.

naive_model_buffered_lynx <- glm(coyotedetection ~ latitude+UrbanIndustrial+UrbanResidence+RoadHardSurface,
                            family = binomial, data = naive_data_buffered)

##Potential Climate Variables: MCMT, MAT, FFP
##Human Footprint: RoadHardSurface,

# camtrapR occupancy modeling workflow

##Camera operation matrix
camera_operation_dates_filtered<-camera_operation_dates_v2%>%
  filter(!is.na(start_date_time))

camera_operation_dates_filtered<-camera_operation_dates_filtered%>%
  mutate(start_date_time = as.character(start_date_time))%>%
  mutate(end_date_time = as.character(end_date_time))

##Start and end dates for deployments per year
camera_deployment_dates<-camera_operation_dates_filtered %>%
  group_by(location, year)%>%
  summarize(
    start_date_time = min(start_date_time),
    end_date_time = max (end_date_time),
    .groups = "drop"
  )

##Idenitfy gaps or periods of inactivity within a camera's deployment "problems"
problem_matrix <-camera_operation_dates_filtered%>%
  arrange(location, year, start_date_time)%>%
  group_by(location, year)%>%
  mutate(
    Problem1_from = lag(end_date_time) + 1,
    Problem1_to = start_date_time -1
  ) %>%
  filter(!is.na(Problem1_from)& Problem1_from <= Problem1_to)%>%
  select(location, year, Problem1_from, Problem1_to)

##Combine problems with deployment data
camera_deployment_problems<-camera_deployment_dates%>%
  left_join(problem_matrix, by = c("location", "year"))%>%
              arrange(location, year, start_date_time, Problem1_from)

##Convert dates to character to be compatible with camptrapR
camera_deployment_problems<-camera_deployment_problems%>%
  mutate(start_date_time = as.character(start_date_time))%>%
  mutate(end_date_time = as.character(end_date_time))

##Filter cameras where start and end time are identical
camera_deployment_problems<-camera_deployment_problems%>%
  filter(start_date_time != end_date_time)

#Build a camera operation matrix to fill with detections
camera_operation_matrix<-cameraOperation(camera_deployment_problems,
                                         stationCol = "location",
                                         sessionCol = "year",
                                         setupCol = "start_date_time",
                                         retrievalCol = "end_date_time",
                                         writecsv = FALSE,
                                         hasProblems = TRUE,
                                         dateFormat = "ymd_HMS"
                                         )
##Camera detection history
datetime_check<-as.POSIXct(study_main_reports$image_date_time, format = ymd_HMS)
invalid_entries<-study_main_reports[is.na(datetime_check),]
invalid_entries

###Daylight savings time!! 4 coyote observations fall in the DST hour, so shifting to 3:00 to compensate (I'll figure out a more elegant way to do this later)

##entries 255396,255397,1882366,7326233

study_main_reports[255396,9] <- as.POSIXct("2015-03-08 3:03:29", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
study_main_reports[255397,9]<- as.POSIXct("2015-03-08 3:03:49", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
study_main_reports[1882366,9]<- as.POSIXct("2015-03-08 3:00:01", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
study_main_reports[7326233,9]<- as.POSIXct("2018-03-11 3:15:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

coyotedetectionhistory<-detectionHistory(recordTable = study_main_reports,
                                         camOp = camera_operation_matrix,
                                         stationCol = "location",
                                         speciesCol = "species_common_name",
                                         recordDateTimeCol = "image_date_time",
                                         recordDateTimeFormat = "ymd_HMS",
                                         species = "Coyote",
                                         occasionLength = 1,
                                         day1 = "station",
                                         datesAsOccasionNames = FALSE,
                                         includeEffort = TRUE,
                                         scaleEffort = FALSE,
                                         timeZone = "America/Edmonton",
                                         unmarkedMultFrameInput = TRUE
                                         )
coyotedetections<-coyotedetectionhistory$detection_history

##Survey report
survey_report <- surveyReport(recordTable = study_main_reports,
                              CTtable = camera_operation_dates,
                              camOp = camera_operation_matrix,
                              speciesCol = "species_common_name",
                              stationCol = "location",
                              setupCol = "start_date_time",
                              retrievalCol = "end_date_time",
                              CTDateFormat = "ymd_HMS",
                              recordDateTimeCol = "image_date_time",
                              recordDateTimeFormat = "ymd_HMS"
                              )

##Site Covariates
site_covariates<-study_main_reports%>%
  distinct(location, latitude, longitude)

###Climate


###Human Footprint

###Convert to unmarkedFrameOccu object

site_covariates <- site_covariates[match(rownames(coyotedetections), site_covariates$location), ]

umf <- unmarkedFrameOccu(
  y = coyotedetectionhistory$detection_history,
  siteCovs = sitecovariates
)

#Exploration of overall detections across latitudes

# Filter to carnivores
carnivores <- all_main_reports %>%
  dplyr::filter(species_common_name %in% c("Bear", "Black Bear", "Bobcat", "Canada Lynx",
                                           "Canidae Family", "Common Raccoon", "Cougar",
                                           "Coyote", "Ermine", "Fisher", "Gray Wolf",
                                           "Grizzly Bear", "Marten", "Red Fox", "Striped Skunk",
                                           "Swift fox", "Weasels and Ermine"))

# plot counts
ggplot(carnivores, aes(x = forcats::fct_infreq(species_common_name))) +
  geom_bar() +
  theme_bw() +
  coord_flip()

# plot number of cameras at which species were detected
carnivores %>%
  dplyr::select(species_common_name, location) %>%
  unique() %>%
  dplyr::count(species_common_name) %>%
  dplyr::mutate(percent = n/length(unique(all_main_reports$location))) %>%
  ggplot(aes(x = fct_reorder(species_common_name, percent), y = percent)) +
  geom_col() +
  theme_bw() +
  coord_flip() +
  xlab("Percent of Locations Detected")

# create column for year
carnivores$year <- lubridate::year(carnivores$image_date_time)

# filter out 2014 and 2020
carnivores <- carnivores %>%
  filter(between(year, 2015, 2020))

# mapping detections with camtrapR

maptest<-detectionMaps(CTtable = carnivores,
                       recordTable = carnivores,
                       Xcol= "longitude",
                       Ycol = "latitude",
                       stationCol = "location",
                       speciesCol = "species_common_name",
                       writePNG = FALSE,
                       plotR = TRUE,
                       printLabels = TRUE,
                       richnessPlot = TRUE,
                       addLegend = TRUE)

# coyotes over the years
carnivores %>%
  dplyr::filter(species_common_name == "Coyote") %>%
  dplyr::count(year) %>%
  ggplot(aes(x = year, y = n)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("Number of Coyote Detections")

# coyotes over latitude over the years
count(filter(carnivores, species_common_name == "Coyote"))

carnivores %>%
  dplyr::filter(species_common_name == "Coyote") %>%
  ggplot(aes(x = as.factor(year), y = latitude)) +
  geom_boxplot() +
  geom_point() +
  theme_bw()

#does mean latitude of detection vary?
latitudemodel_coyote <-lm(latitude~year, data = carnivores, subset=(species_common_name == "Coyote"))
summary(latitudemodel_coyote)
latitude_aov <- aov(latitude~as.factor(year), data = carnivores)
summary(latitude_aov)
TukeyHSD(latitude_aov, conf.level = .95)

#marten detections across latitude by year
count(filter(carnivores, species_common_name == "Marten"))

carnivores%>%
  dplyr::filter(species_common_name == "Marten")%>%
  ggplot(aes(x=as.factor(year), y=latitude))+
  geom_boxplot()+
  geom_point()+
  theme_bw()+
  xlab("Number of Marten Detections")

#marten detections in 2017 suddenly across al latitudes, after being limited to higher latitudes. What could be going on?

#fisher detections across latitude by year
count(filter(carnivores, species_common_name == "Fisher"))

carnivores%>%
  dplyr::filter(species_common_name == "Fisher")%>%
  ggplot(aes(x=as.factor(year), y=latitude))+
  geom_boxplot()+
  geom_point()+
  theme_bw()+
  xlab("Number of Fisher Detections")


#red fox detections across latitude by year
carnivores%>%
  dplyr::filter(species_common_name == "Red Fox")%>%
  ggplot(aes(x=as.factor(year), y=latitude))+
  geom_boxplot()+
  geom_point()+
  theme_bw()+
  xlab("Number of Red Fox Detections")

#canada lynx detections across latitude by year
count(filter(carnivores, species_common_name == "Canada Lynx"))

carnivores%>%
  dplyr::filter(species_common_name == "Canada Lynx")%>%
  ggplot(aes(x=as.factor(year), y=latitude))+
  geom_boxplot()+
  geom_point()+
  theme_bw()+
  xlab("Number of Canada Lynx Detections")

#all 3 taxa
carnivores%>%
  dplyr::filter(species_common_name == "Marten" |
                  species_common_name == "Coyote" |
                  species_common_name == "Fisher")%>%
  ggplot(aes(x=as.factor(year), y=latitude, fill = species_common_name))+
  geom_boxplot()+
  geom_smooth(method = "lm")+
  theme_bw()+
  xlab("Number of Detections")

# something is up with 2020, not much data yet?


# where is the meta-data on camera operation dates?

#Species activity
activityDensity(recordTable = carnivores,
                speciesCol = "species_common_name",
                species = "Coyote",
                recordDateTimeCol = "image_date_time",
                recordDateTimeFormat = "ymd HMS",
                PlotR = TRUE,
                writePNG = FALSE)

##Lynx Coyote Overlap
activityOverlap(recordTable = carnivores,
                speciesCol = "species_common_name",
                speciesA = "Canada Lynx",
                speciesB = "Coyote",
                recordDateTimeCol = "image_date_time",
                recordDateTimeFormat = "ymd HMS",
                writePNG = FALSE,
                plotR = TRUE,
                createDir = FALSE,
                pngMaxPix = 1000,
                linecol = c("black", "blue"),
                linewidth = c(5,3),
                linetype = c(1,2),
                olapool = "darkgrey",
                extend = "lightgrey",
                ylim = c(0, 0.25),
                main = paste("Activity overlap between Lynx and Coyote"))

##Marten Coyote Overlap
activityOverlap(recordTable = carnivores,
                speciesCol = "species_common_name",
                speciesA = "Marten",
                speciesB = "Coyote",
                recordDateTimeCol = "image_date_time",
                recordDateTimeFormat = "ymd HMS",
                writePNG = FALSE,
                plotR = TRUE,
                createDir = FALSE,
                pngMaxPix = 1000,
                linecol = c("black", "blue"),
                linewidth = c(5,3),
                linetype = c(1,2),
                olapool = "darkgrey",
                extend = "lightgrey",
                ylim = c(0, 0.25),
                main = paste("Activity overlap between Marten and Coyote"))

##Fisher Coyote Overlap
activityOverlap(recordTable = carnivores,
                speciesCol = "species_common_name",
                speciesA = "Fisher",
                speciesB = "Coyote",
                recordDateTimeCol = "image_date_time",
                recordDateTimeFormat = "ymd HMS",
                writePNG = FALSE,
                plotR = TRUE,
                createDir = FALSE,
                pngMaxPix = 1000,
                linecol = c("black", "blue"),
                linewidth = c(5,3),
                linetype = c(1,2),
                olapool = "darkgrey",
                extend = "lightgrey",
                ylim = c(0, 0.25),
                main = paste("Activity overlap between Fisher and Coyote"))

##red Fox Coyote Overlap
activityOverlap(recordTable = carnivores,
                speciesCol = "species_common_name",
                speciesA = "Red Fox",
                speciesB = "Coyote",
                recordDateTimeCol = "image_date_time",
                recordDateTimeFormat = "ymd HMS",
                writePNG = FALSE,
                plotR = TRUE,
                createDir = FALSE,
                pngMaxPix = 1000,
                linecol = c("black", "blue"),
                linewidth = c(5,3),
                linetype = c(1,2),
                olapool = "darkgrey",
                extend = "lightgrey",
                ylim = c(0, 0.25),
                main = paste("Activity overlap between Red Fox and Coyote"))

##WolfCoyote Overlap
activityOverlap(recordTable = carnivores,
                speciesCol = "species_common_name",
                speciesA = "Gray Wolf",
                speciesB = "Coyote",
                recordDateTimeCol = "image_date_time",
                recordDateTimeFormat = "ymd HMS",
                writePNG = FALSE,
                plotR = TRUE,
                createDir = FALSE,
                pngMaxPix = 1000,
                linecol = c("black", "blue"),
                linewidth = c(5,3),
                linetype = c(1,2),
                olapool = "darkgrey",
                extend = "lightgrey",
                ylim = c(0, 0.25),
                main = paste("Activity overlap between Gray Wolf and Coyote"))

##Now the real question, can we partition mesocarnivore to see if diel activity varies in areas with or without coyotes?

##CapTrapR multi-species
study_carnivores <- study_main_reports %>%
  dplyr::filter(species_common_name %in% c("Canada Lynx", "Coyote", "Fisher", "Gray Wolf", "Marten"))

DetHist_list <- lapply(unique(study_carnivores$species_common_name), FUN = function(x) {
  detectionHistory(
    recordTable         = study_carnivores,
    camOp                = camera_operation_matrix,
    stationCol           = "location",
    speciesCol           = "species_common_name",
    recordDateTimeCol    = "image_date_time",
    species              = x,     # this gets modifies by lapply
    occasionLength       = 1,
    day1                 = "station",
    datesAsOccasionNames = FALSE,
    includeEffort        = TRUE,
    scaleEffort          = FALSE,
    timeZone             = "America/Edmonton",
    unmarkedMultFrameInput = TRUE
  )}
)

names(DetHist_list) <- unique(study_carnivores$species_common_name)

ylist <- lapply(DetHist_list, FUN = function(x) x$detection_history)

sitecovariates$location<-as.factor(sitecovariates$location)

data_list <- list(ylist = ylist,
                  siteCovs = sitecovariates,
                  obsCovs = list(effort = DetHist_list[[1]]$effort))

modelfile1 <- tempfile(fileext = ".txt")

data_list$siteCovs <- na.omit(data_list$siteCovs)

mod.jags <- communityModel(data_list,
                           occuCovs = list(fixed = "latitude", ranef = "location"),
                           detCovsObservation = list(fixed = "effort"),
                           intercepts = list(det = "ranef", occu = "ranef"),
                           modelFile = modelfile1)
