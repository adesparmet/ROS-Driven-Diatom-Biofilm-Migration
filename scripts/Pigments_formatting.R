#Library packages ####
#install.packages("readr")
library(readr)
#install.packages("stringr")
library(stringr)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("here")
library(here)

# Data import ####

path <- here("data", "Pigments_data", "Chromato_integration")

files <- list.files(path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE) # .txt file list
data_list <- list() # Create new list

# Extract data
for (file in files) {
  data <- read_delim(file, delim = ";", skip = 2, col_types = cols(.default = col_character()))
  
  data <- data %>%
    mutate(Peak = as.numeric(Peak), Area = as.numeric(Area)) %>%
    select(Peak, Area) # Numeric data
  
  data <- data %>%
    mutate(File = basename(file) %>% tools::file_path_sans_ext()) # File name
  
  data_list[[file]] <- data
}

df_combined_pig <- bind_rows(data_list) # Combined into one df

# Pigment identifiChlation ####

#412nm
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^412", File), "Chla-l1", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 2 & grepl("^412", File), "Chla-l3", Peak))

#431nm
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^431", File), "Chla-l2", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 2 & grepl("^431", File), "Chla-allomer", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 3 & grepl("^431", File), "Chla", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 4 & grepl("^431", File), "Chla-epimer", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 5 & grepl("^431", File), "Chla-l4", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 6 & grepl("^431", File), "Chla-l5", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 7 & grepl("^431", File), "Chla-l6", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 8 & grepl("^431", File), "Chla-l7", Peak))

#441
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^441", File), "N", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 2 & grepl("^441", File), "V", Peak))

#446
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^446", File), "Chlc2", Peak))

#448
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^448", File), "Fl1", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 2 & grepl("^448", File), "F", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 3 & grepl("^448", File), "Fl2", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 4 & grepl("^448", File), "Fl3", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 5 & grepl("^448", File), "Fl4", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 6 & grepl("^448", File), "Fl5", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 7 & grepl("^448", File), "Fl6", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 8 & grepl("^448", File), "Dd", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 9 & grepl("^448", File), "Fl7", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 10 & grepl("^448", File), "A", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 11 & grepl("^448", File), "Uc1", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 12 & grepl("^448", File), "Uc2", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 13 & grepl("^448", File), "Uc3", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 14 & grepl("^448", File), "Uc4", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 15 & grepl("^448", File), "L", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 16 & grepl("^448", File), "Z", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 17 & grepl("^448", File), "Uc5", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 18 & grepl("^448", File), "Uc6", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 19 & grepl("^448", File), "Uc7", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 20 & grepl("^448", File), "Uc8", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 21 & grepl("^448", File), "Uc9", Peak))

#454
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^454", File), "Dt", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 2 & grepl("^454", File), "βε", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 3 & grepl("^454", File), "ββ", Peak))

#663
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 1 & grepl("^663", File), "Pda1", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 2 & grepl("^663", File), "Pda2", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 3 & grepl("^663", File), "Pda3", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 4 & grepl("^663", File), "Phal1", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 5 & grepl("^663", File), "Pha", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 6 & grepl("^663", File), "Phal2", Peak))
df_combined_pig <- df_combined_pig %>%
  mutate(Peak = ifelse(Peak == 7 & grepl("^663", File), "Pya", Peak))
# Coefficient standard pigment director ####

# Chla-libration curve of DHI standard pigments (Suppl MM)
coeff_Chla412<-6809.6001 #Chla at 412nm
coeff_Chla431<-8393.2274 #Chla at 431nm
coeff_Chla663<-7043.3959 #Chla at 663nm
coeff_L<-17820.4549 #Lutein at 448nm
coeff_Chlc2<- 28129.7705 #Chlc2 at 446nm
coeff_V<- 13113.8819 #Violaxanthin at 441nm
coeff_F<-12985.4446 #Fucoxanthin at 448nm
coeff_B<-19562.5542 #Beta-Chlarotene at 454nm
coeff_Dt<-21162.9569 #Diatoxanthin at 454nm
coeff_Dd<-17697.5559 #Diadinoxanthin at 448nm

# Area-to-weight pigment conversion ####

df_combined_pig <- df_combined_pig %>%mutate(µg = NA) #Creation of µg column

df_combined_pig[df_combined_pig$Peak == "Chlc2", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chlc2", "Area"]/coeff_Chlc2
df_combined_pig[df_combined_pig$Peak == "F", "µg"] <- df_combined_pig[df_combined_pig$Peak == "F", "Area"]/coeff_F
df_combined_pig[df_combined_pig$Peak == "N", "µg"] <- df_combined_pig[df_combined_pig$Peak == "N", "Area"]/ coeff_V
df_combined_pig[df_combined_pig$Peak == "V", "µg"] <- df_combined_pig[df_combined_pig$Peak == "V", "Area"]/ coeff_V
df_combined_pig[df_combined_pig$Peak == "Dd", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Dd", "Area"]/ coeff_Dd
df_combined_pig[df_combined_pig$Peak == "A", "µg"] <- df_combined_pig[df_combined_pig$Peak == "A", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Dt", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Dt", "Area"]/ coeff_Dt
df_combined_pig[df_combined_pig$Peak == "L", "µg"] <- df_combined_pig[df_combined_pig$Peak == "L", "Area"]/ coeff_L
df_combined_pig[df_combined_pig$Peak == "Z", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Z", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Chla", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla", "Area"]/ coeff_Chla431
df_combined_pig[df_combined_pig$Peak == "Pha", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Pha", "Area"]/ coeff_Chla663
df_combined_pig[df_combined_pig$Peak == "βε", "µg"] <- df_combined_pig[df_combined_pig$Peak == "βε", "Area"]/ coeff_B
df_combined_pig[df_combined_pig$Peak == "ββ", "µg"] <- df_combined_pig[df_combined_pig$Peak == "ββ", "Area"]/ coeff_B
df_combined_pig[df_combined_pig$Peak == "Pya", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Pya", "Area"]/ coeff_Chla663

df_combined_pig[df_combined_pig$Peak == "Fl1", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl1", "Area"]/coeff_F
df_combined_pig[df_combined_pig$Peak == "Fl2", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl2", "Area"]/coeff_F
df_combined_pig[df_combined_pig$Peak == "Fl3", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl3", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Fl4", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl4", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Fl5", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl5", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Fl6", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl6", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Fl7", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Fl7", "Area"]/ coeff_F

df_combined_pig[df_combined_pig$Peak == "Uc1", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc1", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc2", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc2", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc3", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc3", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc4", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc4", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc5", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc5", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc6", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc6", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc7", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc7", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc8", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc8", "Area"]/ coeff_F
df_combined_pig[df_combined_pig$Peak == "Uc9", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Uc9", "Area"]/ coeff_F

df_combined_pig[df_combined_pig$Peak == "Chla-l1", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l1", "Area"]/ coeff_Chla412
df_combined_pig[df_combined_pig$Peak == "Chla-l2", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l2", "Area"]/coeff_Chla431
df_combined_pig[df_combined_pig$Peak == "Chla-l3", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l3", "Area"]/ coeff_Chla412
df_combined_pig[df_combined_pig$Peak == "Chla-l4", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l4", "Area"]/ coeff_Chla431
df_combined_pig[df_combined_pig$Peak == "Chla-l5", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l5", "Area"]/ coeff_Chla431
df_combined_pig[df_combined_pig$Peak == "Chla-l6", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l6", "Area"]/ coeff_Chla431
df_combined_pig[df_combined_pig$Peak == "Chla-l7", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-l7", "Area"]/ coeff_Chla431

df_combined_pig[df_combined_pig$Peak == "Chla-allomer", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-allomer", "Area"]/ coeff_Chla431
df_combined_pig[df_combined_pig$Peak == "Chla-epimer", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Chla-epimer", "Area"]/ coeff_Chla431

df_combined_pig[df_combined_pig$Peak == "Pda1", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Pda1", "Area"]/ coeff_Chla663
df_combined_pig[df_combined_pig$Peak == "Pda2", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Pda2", "Area"]/ coeff_Chla663
df_combined_pig[df_combined_pig$Peak == "Pda3", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Pda3", "Area"]/ coeff_Chla663

df_combined_pig[df_combined_pig$Peak == "Phal1", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Phal1", "Area"]/ coeff_Chla663
df_combined_pig[df_combined_pig$Peak == "Phal2", "µg"] <- df_combined_pig[df_combined_pig$Peak == "Phal2", "Area"]/ coeff_Chla663

# Calculation of pigment concentration per gram of dry sediment (df_pigment and df_pigment_perc) ####

#Calculation of pigment weight in the extraction volume (2ml methanol)
df_combined_pig <- df_combined_pig %>%mutate(w2ml = NA)
df_combined_pig$w2ml <- ((df_combined_pig$µg * 10)*2) # x10 beChlause the injection volume method was 100 µl and multiplied by 2 beChlause the final volume was = 2ml meoh

#Standardization per gram of dry matter
df_combined_pig <- df_combined_pig %>%mutate(µg.sed = NA)

VLL_WSB_1<-	0.0103 #Weighed dry matter (gr)
VLL_WSB_2<-	0.0118
VLL_WSB_3<-	0.0117
VLL_WSB_4<-	0.0113
VLL_WSB_5<-	0.0112
HL_WSB_1<-	0.0111
HL_WSB_2<-	0.0116
HL_WSB_3<-	0.0108
HL_WSB_4<-	0.0114
HL_WSB_5<-	0.0112
PL_WSB_1<-	0.0109
PL_WSB_2<-	0.0116
PL_WSB_3<-	0.0115
PL_WSB_4<-	0.0113
PL_WSB_5<-	0.0109
HP_WSB_1<-	0.0111
HP_WSB_2<-	0.0113
HP_WSB_3<-	0.0114
HP_WSB_4<-	0.0105
HP_WSB_5<-	0.0112
VLL_SB_1<-	0.0406
VLL_SB_2<-	0.0414
VLL_SB_3<-	0.0417
VLL_SB_4<-	0.0399
VLL_SB_5<-	0.0398
HL_SB_1<-	0.0417
HL_SB_2<-	0.0410
HL_SB_3<-	0.0407
HL_SB_4<-	0.0396
HL_SB_5<-	0.0400
PL_SB_1<-	0.0417
PL_SB_2<-	0.0419
PL_SB_3<-	0.0405
PL_SB_4<-	0.0405
PL_SB_5<-	0.0408
HP_SB_1<-	0.0407
HP_SB_2<-	0.0410
HP_SB_3<-	0.0408
HP_SB_4<-	0.0403
HP_SB_5<-	0.0418

df_combined_pig <- df_combined_pig %>%
  mutate(µg.sed = case_when(
    grepl("VLL_SB_1", File) ~ w2ml / VLL_SB_1, #Concentrations Calculation
    grepl("VLL_SB_2", File) ~ w2ml / VLL_SB_2,
    grepl("VLL_SB_3", File) ~ w2ml / VLL_SB_3,
    grepl("VLL_SB_4", File) ~ w2ml / VLL_SB_4,
    grepl("VLL_SB_5", File) ~ w2ml / VLL_SB_5,
    grepl("VLL_WSB_1", File) ~ w2ml / VLL_WSB_1,
    grepl("VLL_WSB_2", File) ~ w2ml / VLL_WSB_2,
    grepl("VLL_WSB_3", File) ~ w2ml / VLL_WSB_3,
    grepl("VLL_WSB_4", File) ~ w2ml / VLL_WSB_4,
    grepl("VLL_WSB_5", File) ~ w2ml / VLL_WSB_5,
    grepl("HP_SB_1", File) ~ w2ml / HP_SB_1,
    grepl("HP_SB_2", File) ~ w2ml / HP_SB_2,
    grepl("HP_SB_3", File) ~ w2ml / HP_SB_3,
    grepl("HP_SB_4", File) ~ w2ml / HP_SB_4,
    grepl("HP_SB_5", File) ~ w2ml / HP_SB_5,
    grepl("HP_WSB_1", File) ~ w2ml / HP_WSB_1,
    grepl("HP_WSB_2", File) ~ w2ml / HP_WSB_2,
    grepl("HP_WSB_3", File) ~ w2ml / HP_WSB_3,
    grepl("HP_WSB_4", File) ~ w2ml / HP_WSB_4,
    grepl("HP_WSB_5", File) ~ w2ml / HP_WSB_5,
    grepl("PL_SB_1", File) ~ w2ml / PL_SB_1,
    grepl("PL_SB_2", File) ~ w2ml / PL_SB_2,
    grepl("PL_SB_3", File) ~ w2ml / PL_SB_3,
    grepl("PL_SB_4", File) ~ w2ml / PL_SB_4,
    grepl("PL_SB_5", File) ~ w2ml / PL_SB_5,
    grepl("PL_WSB_1", File) ~ w2ml / PL_WSB_1,
    grepl("PL_WSB_2", File) ~ w2ml / PL_WSB_2,
    grepl("PL_WSB_3", File) ~ w2ml / PL_WSB_3,
    grepl("PL_WSB_4", File) ~ w2ml / PL_WSB_4,
    grepl("PL_WSB_5", File) ~ w2ml / PL_WSB_5,
    grepl("HL_SB_1", File) ~ w2ml / HL_SB_1,
    grepl("HL_SB_2", File) ~ w2ml / HL_SB_2,
    grepl("HL_SB_3", File) ~ w2ml / HL_SB_3,
    grepl("HL_SB_4", File) ~ w2ml / HL_SB_4,
    grepl("HL_SB_5", File) ~ w2ml / HL_SB_5,
    grepl("HL_WSB_1", File) ~ w2ml / HL_WSB_1,
    grepl("HL_WSB_2", File) ~ w2ml / HL_WSB_2,
    grepl("HL_WSB_3", File) ~ w2ml / HL_WSB_3,
    grepl("HL_WSB_4", File) ~ w2ml / HL_WSB_4,
    grepl("HL_WSB_5", File) ~ w2ml / HL_WSB_5,
    TRUE ~ NA_real_
  ))

# Matrix of pigment concentration per gram of dry sediment
df_pigment<- df_combined_pig %>%
  mutate(File = str_replace(File, ".*txt_", "")) %>% 
  mutate(File = str_remove(File, "_0524$")) %>%
  select(File, Peak, µg.sed) %>%
  pivot_wider(
    names_from = Peak,
    values_from = µg.sed, 
    values_fill = list(µg.sed = NA))

# Add metadata
df_pigment <- df_pigment %>%mutate(Treatment = str_extract(File, "^[^_]+_[^_]+"))
df_pigment <- df_pigment %>%mutate(Conditions = str_extract(File, "[^_]+"))
df_pigment <- df_pigment %>%mutate(Groups = str_extract(File, "(?<=_)[^_]+"))
df_pigment <- df_pigment %>%select(File, Treatment, Conditions, Groups, everything())

# Percentage conversion
num_data <- df_pigment %>%select(-File, -Treatment, -Conditions, -Groups)
df_pigment_perc <- df_pigment %>%mutate(across(where(is.numeric), ~ . / rowSums(num_data)) * 100)

# Grouped df_pigment (gdf_pigment and gdf_pigment_perc) ####

# Chla_derivative group
gdf_pigment <- df_pigment %>%
  mutate(`Chla-derivatives`= rowSums(select(.,starts_with("Chla-l"), starts_with("Chla-a"),starts_with("Chla-e")))) %>%
  select(-starts_with("Chla-l"), -starts_with("Chla-a"), -starts_with("Chla-e")) #Chla like + allo + epi

# Pheopigment group
gdf_pigment <- gdf_pigment %>%
  mutate(`Pheopigments`= rowSums(select(.,starts_with("Pda"), starts_with("Pha"), starts_with("Pya")))) %>%
  select(-starts_with("Pda"), -starts_with("Pha"), -starts_with("Pya"))

# Transformation in percentage
cnum_data <- gdf_pigment %>%select(-File, -Treatment, -Conditions, -Groups)
gdf_pigment_perc <- gdf_pigment %>%mutate(across(where(is.numeric), ~ . / rowSums(cnum_data)) * 100)

