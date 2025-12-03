#Library ####
#install.packages("jsonlite")
library(jsonlite)
#install.packages("lubridate")
library(lubridate)
#install.packages("hms")
library(hms)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
# Stress migration data ####

File_HL <- "1_Data_formatting/data/Photophysiology_data/Stress_migration_behavior/20240219_HL.txt"

File_HP <- "1_Data_formatting/data/Photophysiology_data/Stress_migration_behavior/20240514_HP.txt"

File_PL <- "1_Data_formatting/data/Photophysiology_data/Stress_migration_behavior/20240216_PL.txt"


# Txt to JSON
File_HL <- readLines(File_HL)
File_HL_data <- fromJSON(paste(File_HL, collapse = ""))

File_HP <- readLines(File_HP)
File_HP_data <- fromJSON(paste(File_HP, collapse = ""))

File_PL <- readLines(File_PL)
File_PL_data <- fromJSON(paste(File_PL, collapse = ""))

# `F0`
File_HL_data$`F0` <- File_HL_data$flash - File_HL_data$backgr
File_HP_data$`F0` <- File_HP_data$flash - File_HP_data$backgr
File_PL_data$`F0` <- File_PL_data$flash - File_PL_data$backgr

# Add stress HL and DA metadata
n_hl1 <- 1:48
n_hl2 <- 49:96
n_hl3 <- 97:248
stress_pattern_hl1 <- rep(c(rep("DA", 8), rep("HL", 8)), length.out = length(n_hl1))
stress_pattern_hl2 <- rep("DA", length(n_hl2)) # No HL measurment during HL stress application
stress_pattern_hl3 <- rep(c(rep("HL", 8), rep("DA", 8)), length.out = length(n_hl3)) # 8 mesures
stress_pattern_hl <- c(stress_pattern_hl1, stress_pattern_hl2, stress_pattern_hl3)
File_HL_data$Conditions<- stress_pattern_hl

# Add stress HP and DA metadata
n_hp <- nrow(File_HP_data)
stress_pattern_hp <- rep(c(rep("HP", 5), rep("DA", 5)), length.out = n_hp) # 5 mesures
File_HP_data$Conditions <- stress_pattern_hp

# Add stress PL and DA metadata
n_pl <- nrow(File_PL_data)
stress_pattern_pl <- rep(c(rep("DA", 8), rep("PL", 8)), length.out = n_pl)  # 8 mesures
File_PL_data$Conditions <- stress_pattern_pl

# Time format
File_HL_data$time <- sub(" .*", "", File_HL_data$time)
File_HP_data$time <- sub(" .*", "", File_HP_data$time)
File_PL_data$time <- sub(" .*", "", File_PL_data$time)

File_HL_data$time <- hms::as_hms(File_HL_data$time)
File_HP_data$time <- hms::as_hms(File_HP_data$time)
File_PL_data$time <- hms::as_hms(File_PL_data$time)

# Bug time correction for HL
vl_72 <- File_HL_data$time[72]
File_HL_data$time[73:248] <- File_HL_data$time[73:248] + as.duration(vl_72)

# Time adjustement
vl_1_hl <- File_HL_data$time[1]
File_HL_data$time <- as_hms(File_HL_data$time - as.duration(vl_1_hl))

vl_1_hp <- File_HP_data$time[1]
File_HP_data$time <- as_hms(File_HP_data$time - as.duration(vl_1_hp))

vl_1_pl <- File_PL_data$time[1]
File_PL_data$time <- as_hms(File_PL_data$time - as.duration(vl_1_pl))

# Minute conversions
File_HL_data <- File_HL_data %>%
  mutate(minutes = as.numeric(time) / 60)
File_HP_data <- File_HP_data %>%
  mutate(minutes = as.numeric(time) / 60)
File_PL_data <- File_PL_data %>%
  mutate(minutes = as.numeric(time) / 60)