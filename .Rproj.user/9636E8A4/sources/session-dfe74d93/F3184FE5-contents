#Library packages ####
#install.packages("jsonlite")
library(jsonlite)
#install.packages("utils")
library(utils)
#install.packages("devtools")
library(devtools)
#install_github("bmjesus/thejesus")
library(thejesus)
#install.packages("tidyr")
library(tidyr)
#install.packages("dplyr")
library(dplyr)
#install.packages("gt")
library(gt)
#install.packages("here")
library(here)

# Import RLC data ####
txt_content_HL <- readLines(here("data", "Photophysiology_data", "Photosynthetic_parameters", "20240501_HL.txt"), warn = FALSE)

txt_content_PL <- readLines(here("data", "Photophysiology_data", "Photosynthetic_parameters", "20240502_PL.txt"), warn = FALSE)

txt_content_HP <- readLines(here("data", "Photophysiology_data", "Photosynthetic_parameters", "20240503_HP.txt"), warn = FALSE)

# Merge raw data files
txt_content_HL <- txt_content_HL[txt_content_HL != "]"]
txt_content_PL <- txt_content_PL[txt_content_PL != "["]
txt_content_HL[length(txt_content_HL)] <- sub("\\s*}\\s*$", "},", txt_content_HL[length(txt_content_HL)])
txt_content_PL <- txt_content_PL[txt_content_PL != "]"]
txt_content_HP <- txt_content_HP[txt_content_HP != "["]
txt_content_PL[length(txt_content_PL)] <- sub("\\s*}\\s*$", "},", txt_content_PL[length(txt_content_PL)])
txt_content <- c(txt_content_HL, txt_content_PL, txt_content_HP)

# Data extraction
json_data <- paste0(txt_content, collapse = "") 

df <- as.data.frame(t(fromJSON(json_data)[, c("fo", "fm", "fm_l1", "fm_l2", "fm_l3", "fm_l4", "fm_l5", "fm_l6", "fm_l7", "ft_l1", "ft_l2", "ft_l3", "ft_l4", "ft_l5", "ft_l6", "ft_l7", "qy_max", "qy_l1", "qy_l2", "qy_l3", "qy_l4", "qy_l5", "qy_l6", "qy_l7")])) # df creation
df_transposed <- t(df)
df_transposed <- as.data.frame(df_transposed)

colonnes_numeriques <- sapply(df_transposed, is.numeric)
df_transposed[, colonnes_numeriques] <- lapply(df_transposed[, colonnes_numeriques], function(x) gsub(",", ".", x))

# Add metadata
metadata_photophy <- read.csv(here("data", "Photophysiology_data", "Photosynthetic_parameters", "metadata_photophy.csv"), header = TRUE, sep = ";") # metadata file

# Merge metadata and df_transposed
lc.data <- merge(df_transposed, metadata_photophy, by = "row.names", all = TRUE)
lc.data <- lc.data[, -1]
lc.data[, 1:24] <- lapply(lc.data[, 1:24], as.numeric)

# rETR calculation ####

PSII<-data.frame(`Fq_Fm`=(lc.data[,2]-lc.data[,1])/lc.data[,2])

for( i in 1:7){
  PSII[,i+1]<-(lc.data[,i+2]-lc.data[,i+9])/lc.data[,i+2] 
} #PSII efficiency calculation

lights<-c(10,20,50,100,300,500,1000) # RLC program light steps

rETR<-t(t(PSII)*c(0,lights)*0.5) # rETR calculation, *0.5 assuming 50/50 photons distribution between PSI/PSII
colnames(rETR)<-c("rETR 0",paste("rETR",1:7))
rETR<-data.frame(rETR) # df creation

# Retrieve rETRm_obs
rETR$`rETRm obs` <- apply(rETR[, 1:8], 1, max, na.rm = TRUE) # Retrieve rETRm observed
all(apply(rETR[, c(8,9)], 1, function(x) x[1] == x[2])) # All rETRm obs at 1000 PPFD?

# Photosynthetic parameters calculation (Eilers&Peeters) ####

source(here("data", "Photophysiology_data", "modele_non_lineaire_Eilers&Peeters.R")) #Eilers&Peeters script source

# Model calculation
res <- list()
for(i in 1:nrow(rETR)){  
  res[[i]] = fit_ep(light = c(0, lights), as.numeric(rETR[i, 1:8]))  
}

df_Pparameters<-do.call("rbind", res) # df creation
df_Pparameters$F0<-lc.data$fo # Add F0 values
df_Pparameters$"Fq'/Fm' VLL"<-lc.data$qy_max # Add Fq'/Fm' values
df_Pparameters$Name<-lc.data$Name # Add metadata
df_Pparameters$No<-lc.data$No
df_Pparameters$Time<-lc.data$Time
df_Pparameters$Treatment<-lc.data$Treatment

# Rename
df_Pparameters <- df_Pparameters %>%rename(rETRm_model = rETRm)
df_Pparameters <- df_Pparameters %>%rename("F0" = F0)

df_Pparameters <- df_Pparameters %>%
  mutate(`F0` = ifelse(grepl("^WSB_", Treatment), NA, `F0`)) #Useless values because no migration

# Add "Day" variable metadata
df_Pparameters$Day <- ifelse(grepl("_HL_|_VLL_1|_VLL_2|_VLL_3", df_Pparameters$Name), 1, 
                      ifelse(grepl("_PL_|_VLL_4|_VLL_5|_VLL_6", df_Pparameters$Name), 2, 
                     ifelse(grepl("_HP_|_VLL_7|_VLL_8|_VLL_9", df_Pparameters$Name), 3, NA)))

# NPQ, YNPQ, YNO calculations ####

# NPQ (1000)
for (i in 1:7) {
  fm_prime <- lc.data[[paste0("fm_l", i)]]
  fm_dark  <- lc.data$fm
  lc.data[[paste0("NPQ_l", i)]] <- (fm_dark - fm_prime) / fm_prime
}

# YNO (max_obs)
for (i in 1:7) {
  f0_prime <- lc.data[[paste0("ft_l", i)]]  
  fm_dark  <- lc.data$fm                    
  lc.data[[paste0("YNO_l", i)]] <- f0_prime / fm_dark
}


# YNPQ (1000)
for (i in 1:7) {
  fm_prime <- lc.data[[paste0("fm_l", i)]]  
  f0_prime <- lc.data[[paste0("ft_l", i)]]  
  fm_dark  <- lc.data$fm                     
  
  lc.data[[paste0("YNPQ_l", i)]] <- (f0_prime / fm_prime) - (f0_prime / fm_dark)
}


# YNPQ + YNO + QY = 1 ?
for (i in 1:7) {
  fq_prime <- lc.data[[paste0("qy_l", i)]]  
  yno_prime <- lc.data[[paste0("YNO_l", i)]]  
  ynpq_prime <- lc.data[[paste0("YNPQ_l", i)]]                      
  
  lc.data[[paste0("eff_tot", i)]] <- (fq_prime + yno_prime + ynpq_prime)
}

# Quality check
check_efftot <- function(df) {
  cols_to_check <- df[, 50:56]
  all(cols_to_check >= 0.99 & cols_to_check <= 1.01, na.rm = TRUE)
}
check_efftot(lc.data) # All total values between 0.99 and 1.01?

# Clean negative artefacts
lc.data <- lc.data %>% mutate(across(
  .cols = matches("^(NPQ|YNPQ|YNO)"),
  .fns = ~ ifelse(.x < 0, 0, .x)))

# Retrieve npq, ynpq and yno max

#NPQ
lc.data$NPQ1000 <- apply(lc.data[, 29:35], 1, max, na.rm = TRUE) # Retrieve NPQm obs
all(apply(lc.data[, c(35,57)], 1, function(x) x[1] == x[2])) # All NPQ max obs at 1000 PPFD?

#YNPQ
lc.data$YNPQ1000 <- apply(lc.data[, 43:49], 1, max, na.rm = TRUE) # Retrieve YNPQm obs
all(apply(lc.data[, c(49,58)], 1, function(x) x[1] == x[2])) # All YNPQ max obs at 1000 PPFD?

#YNO
lc.data$`YNOm obs` <- apply(lc.data[, 36:42], 1, max, na.rm = TRUE) # Retrieve YNOm obs
all(apply(lc.data[, c(42,59)], 1, function(x) x[1] == x[2])) # All YNO max obs at 1000 PPFD?

# Merge df photosynthetic parameters (df_Pparameters) ####

lc.data$`rETRm obs` <- rETR$`rETRm obs`
df_Pparameters <- merge(lc.data, df_Pparameters, by = c("Name", "Time", "Treatment", "No"), all = TRUE)

# Check merged
all(apply(df_Pparameters[, c(21,67)], 1, function(x) x[1] == x[2])) # same fq/fm check merged?

# Df formatting
df_Pparameters <- df_Pparameters %>% select(Name, Time, Treatment, No, Day, everything())