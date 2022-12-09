###
###
#' 
#' Script for:
#' TITLE
#' McGlade  et al.
#' Preprint: 
#' 
#' Latest update: 
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Preparing dataset for incubation analysis
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, 
               lubridate, dplyr, tidyr,
               lme4, performance, rptR,
               suncalc,
               ggplot2, extrafont)
loadfonts()


##
##### functions #####
##
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

##
##
##### data #####
##
##
data <- read.csv("./data/light_exp_dataset.csv")
head(data)

##
## new columns
data$inc_start_aprildays <- yday(dmy(data$inc_start)) # incubation start in days after January 1
data$days_to_hatch <- as.numeric((dmy(data$hatching_date) - dmy(data$date))*-1)

##
## clean dataset
data <- data %>% 
  filter(days_to_hatch < 0) %>%  # only one observation for hatching day
filter(days_to_hatch >= -16)

##
##
##### Sunrise and sunset calculation for start of incubation and date of observation #####
##
##
## site coordinates
# table with coordinates for the different sites
coor_df <- data.frame(site = c("KG", "GC", "CASH", "SAL", "SCENE"),
                      lat = c(55.868230, 55.904400, 56.110330, 56.124392, 56.130595),
                      lon = c(-4.282496, -4.322740, -4.576820, -4.601727, -4.614817))


# new columns
data$sunrise_inc_start <- ymd_hms(NA)
data$sunset_inc_start <- ymd_hms(NA)

for(i in 1:nrow(data)){
  location <- data$site[i]
  suntimes  <- getSunlightTimes(date = dmy(data$inc_start)[i], 
                                lat = coor_df$lat[coor_df$site == location],
                                lon = coor_df$lon[coor_df$site == location],
                                keep = c("sunrise", "sunset"), 
                                tz = "Europe/London")
  data$sunrise_inc_start[i] <- suntimes$sunrise
  data$sunset_inc_start[i] <- suntimes$sunset
}

data$sunrise_inc_start_dec <- hour(ymd_hms(data$sunrise_inc_start)) + (minute(ymd_hms(data$sunrise_inc_start))/60)
data$sunset_inc_start_dec <- hour(ymd_hms(data$sunset_inc_start)) + (minute(ymd_hms(data$sunset_inc_start))/60)

##
##
##### Select variables for analysis and save data #####
##
##
data_save <- data %>% 
  select(year, box, date, site, days_to_hatch, inc_start_aprildays, type,
         clutch_size, meantemp,
         first_offbout, 
         last_onbout,
         night_var = night.var,
         activity_onset_relative,
         activity_end_relative)

saveRDS(object = data_save, file = "./data/clean_data_analysis.RDS")





