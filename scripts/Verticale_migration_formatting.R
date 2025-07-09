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
#install.packages("here")
library(here)

# Endogen migration data ####

Engogenous_mig <- here("data", "Photophysiology_data", "Endogen_cycle", "Endogen_migration.txt")

# Txt to JSON
Engogenous_mig <- readLines(Engogenous_mig)
Engogenous_mig <- fromJSON(paste(Engogenous_mig, collapse = ""))

# F0
Engogenous_mig$`F0` <- Engogenous_mig$flash - Engogenous_mig$backgr

# Création du df de travail
Engogenous_mig <- data.frame(`F0` = Engogenous_mig$`F0`, check.names = FALSE)

# Add time
Time <- seq(as.POSIXct("2023-12-22 16:00:00"), by = "30 min", length.out = 609)
data_Time <- data.frame(Time)
Engogenous_mig<-cbind(Engogenous_mig, data_Time["Time"])

# Add Daylight True/False data
hours <- as.numeric(format(data_Time$Time, "%H"))
minutes <- as.numeric(format(data_Time$Time, "%M"))

# Categories
data_Time$Daylight <- (hours > 8 | (hours == 8 & minutes >= 46)) & 
  (hours < 17 | (hours == 16 & minutes == 59)) #set to sunset/sunrise (rise 8:46, set 5:00)
Engogenous_mig<-cbind(Engogenous_mig, data_Time["Daylight"])

# Importation des data de marée
Sealevel_data <- read.table(here("data", "Photophysiology_data", "Endogen_cycle", "Shom_tidal_data.txt"),
                            header = FALSE,
                            sep = "\t")

# df format
Sealevel_data <- Sealevel_data %>%
  separate(V1, into = c("Time", "SeaLevel", "Value"), sep = ";")

Sealevel_data <- subset(Sealevel_data, select = -Value)

# Temporal format
Sealevel_data$Time <- as.POSIXct(Sealevel_data$Time, format = "%d/%m/%Y %H:%M:%S")

# UTC to UTC+1 (Paris) conversion
Sealevel_data$Time <- Sealevel_data$Time + 3600 

# Time data sealevel settings
start_time <- as.POSIXct("2023-12-22 16:00:00", format = "%Y-%m-%d %H:%M:%S")
end_time <- as.POSIXct("2024-01-04 08:00:00", format = "%Y-%m-%d %H:%M:%S")
Sealevel_data <- Sealevel_data[Sealevel_data$Time >= start_time & Sealevel_data$Time <= end_time, ]
Sealevel_data$SeaLevel<-as.numeric(Sealevel_data$SeaLevel)

# Merge
Engogenous_mig_df <- merge(Engogenous_mig, Sealevel_data, by = "Time", all = TRUE)
Engogenous_mig_df <- Engogenous_mig_df[!is.na(Engogenous_mig_df$`F0`), ] # remove na valuesc

# Add hour
Engogenous_mig_df$hour <- format(Engogenous_mig_df$Time, "%H:%M:%S")

# Stress migration data ####

File_HL <- here("data", "Photophysiology_data", "Stress_migration_behavior", "20240219_HL.txt")
File_HP <- here("data", "Photophysiology_data", "Stress_migration_behavior", "20240514_HP.txt")
File_PL <- here("data", "Photophysiology_data", "Stress_migration_behavior", "20240216_PL.txt")

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

# Add stress HL and VLL metadata
n_hl1 <- 1:48
n_hl2 <- 49:96
n_hl3 <- 97:248
stress_pattern_hl1 <- rep(c(rep("VLL", 8), rep("HL", 8)), length.out = length(n_hl1))
stress_pattern_hl2 <- rep("VLL", length(n_hl2)) # No HL measurment during HL stress application
stress_pattern_hl3 <- rep(c(rep("HL", 8), rep("VLL", 8)), length.out = length(n_hl3)) # 8 mesures
stress_pattern_hl <- c(stress_pattern_hl1, stress_pattern_hl2, stress_pattern_hl3)
File_HL_data$Conditions<- stress_pattern_hl

# Add stress HP and VLL metadata
n_hp <- nrow(File_HP_data)
stress_pattern_hp <- rep(c(rep("HP", 5), rep("VLL", 5)), length.out = n_hp) # 5 mesures
File_HP_data$Conditions <- stress_pattern_hp

# Add stress PL and VLL metadata
n_pl <- nrow(File_PL_data)
stress_pattern_pl <- rep(c(rep("VLL", 8), rep("PL", 8)), length.out = n_pl)  # 8 mesures
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