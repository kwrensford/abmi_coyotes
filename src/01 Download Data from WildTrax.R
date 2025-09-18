#-----------------------------------------------------------------------------------------------------------------------

# Title:       Obtain ABMI EH Data for Analysis
# Description: Code for downloading ABMI Ecosystem Health (EH) camera data from the years 2014-2019

# Authors:     Marcus Becker, Kwasi Wrensford, Kaitlyn Gaynor, Melanie Dickie

# Date:        April 2024

#-----------------------------------------------------------------------------------------------------------------------

# Download and attach packages:

# The `wildRtrax` package is used to connect to WildTrax and download camera data
# remotes::install_github('ABbiodiversity/wildRtrax@main')
library(wildRtrax)

# Set up username and password as environment variables
# Note: Has to be 'WT_USERNAME' and 'WT_PASSWORD'
Sys.setenv(WT_USERNAME = "", WT_PASSWORD = "")

# Don't save the credentials in this script and push to Github - not secure.
# I personally use the `keyring` package to store those values for me locally
library(keyring)
Sys.setenv(WT_USERNAME = key_get("WT_USERNAME", keyring = "wildtrax"),
           WT_PASSWORD = key_get("WT_PASSWORD", keyring = "wildtrax"))

# You could also save a small script locally to store these for you
credentials <- "Sys.setenv(WT_USERNAME = '', WT_PASSWORD = '')"
writeLines(credentials, "login.R")
# Then, at the top of your data download script, source the file
source("login.R")

# Authenticate into WildTrax
wt_auth()

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(tidyverse)

# Obtain project ids for ABMI EH 2014-2019 camera projects
eh_projects <- wt_get_download_summary(
  sensor_id = "CAM") |>
  filter(str_detect(project, "Ecosystem Health")) |>
  mutate(year = as.numeric(str_extract_all(project, "\\d+"))) |>
  filter(year < 2020)

# Download data

# Can do one project at a time, e.g., EH 14
eh14_main_report <- wt_download_report(
  # Project ID from dataframe above
  project_id = 205,
  # Specify sensor type (other is 'ARU', which we don't want)
  sensor_id = "CAM",
  # Several report types - "main" gives us mostly all we need
  report = "main"
)

# Or, download all at once into one giant dataframe
# Warning: Could take a little while!

# Vector of project IDs
eh_projects_id <- eh_projects$project_id

all_main_reports <- map_df(
  # Supply all the IDs as a vector
  .x = eh_projects_id,
  # Iterate over them with the download report function
  .f = ~ wt_download_report(
      project_id = .x,
      sensor_id = "CAM",
      report = "main"))

# Download image reports (authentication failure when dowloading image reports, but not for main reports)

wt_get_download_summary(sensor_id = "CAM")

image_report_2014 <- wt_download_report(
  project_id = 205,
  sensor_id = "CAM",
  report = "image_report")

all_image_reports <- map_df(
  .x = eh_projects_id,
  .f = ~ wt_download_report(
    project_id = .x,
    sensor_id = "CAM",
    report = "image_report"))

# Save data for future analysis so we don't have to download each time

# One possibility: save to a `data` folder that is in the .gitignore

write_csv(eh14_main_report, "data/ABMI EH 2014 Main Report.csv")

write_csv(eh14_image_report, "dagta/ABMI EH 2014 Image Report.csv")

# or,
# write_csv(all_main_reports, "data/ABMI EH 2014-2019 Main Reports.csv")

write_csv(all_image_reports, "data/ABMI EH 2014-2019 Image Reports.csv")





