#Library ####
library(readxl)
library(tidyverse)
library(patchwork)
library(dplyr)
library(stringr)
library(plotly)
library(gt)

# Samples total concentration (mg.m-2) ####

# Import pigment concentration data (µg pig/gr dry matter)
df_Qphar <- df_pigment

# Group pigments, derivatives merge
df_Qphar <- df_Qphar %>%
  mutate(`Chlorophyll-a`= rowSums(select(.,starts_with("Chla")))) %>%
  select(-starts_with("Chla"))

df_Qphar <- df_Qphar %>%
  mutate(`Pheophorbide-a`= rowSums(select(.,starts_with("Pda")))) %>%
  select(-starts_with("Pda"))

df_Qphar <- df_Qphar %>%
  mutate(`Pheophytin-a`= rowSums(select(.,starts_with("Pha"), starts_with("Pya")))) %>%
  select(-starts_with("Pha"), -starts_with("Pya"))

df_Qphar <- df_Qphar %>%
  mutate(Fucoxanthin = rowSums(across(c(F, starts_with("Uc"), starts_with("F-l")), as.numeric))) %>%
  select(-F, -starts_with("Uc"), -starts_with("F-l"))

# Conversion de C µg.g -> mg.g-1 
df_Qphar[, 5:18] <- df_Qphar[, 5:18] / 1000 


# Mtot calculation per Petri Dish (mg pig tot/petri dish)

mtot <- "1_Data_formatting/data/Qphar/mtot_samples.xlsx"

df_mtot <- read_excel(mtot, sheet = 1) %>% select(1, 9)

df_Qphar <- df_Qphar %>%
  left_join(df_mtot %>% select(Sample, sum_gr),
            by = c("File" = "Sample")) %>%
  mutate(across(5:18, ~ .x * sum_gr)) %>%
  select(-sum_gr) 


# Formatting 2 df
df_Qphar_Surf <- df_Qphar # Pigment concentration on the Petri Dish surface (mg pig.m-2)
df_Qphar_Vol <- df_Qphar # Pigment concentration in the MPB biovolume (mg pig.m-3)

# MPB (sedimented-diatoms) biovolume in WSB Petri dish 
biofilm_V <- ((pi * (5.5/2)^2)* (0.07332/10))* 1e-6 #Petri Dish diameter = 5.5cm , biofilm height ~ 0.07332 mm (considering a stack of 4 diato which have an average height of 18,33 µm calculated below based on De tommasi et al doi.org/10.1038/s41598-024-56206-y values) 

df_Qphar_Surf[, 5:18] <- (df_Qphar_Surf[, 5:18] / biofilm_V) * 0.00007332 # mg.m-2
df_Qphar_Vol[, 5:18] <- (df_Qphar_Vol[, 5:18] / biofilm_V) # mg.m-3


# Absorption pigment coefficient (m2.mg-1) ####

#doi.org/10.1016/j.dib.2019.103875

SW1 <- "1_Data_formatting/data/Qphar/Clementson and Wojtasiewicz 2019_1.xls"
SW2 <- "1_Data_formatting/data/Qphar/Clementson and Wojtasiewicz 2019_2.xls"
l_SW1 <- read_excel(SW1, sheet = 1)
l_SW2 <- read_excel(SW2, sheet = 1)

l_SW1 <- l_SW1 %>%
  pivot_longer(
    cols = -1,             
    names_to = "Pigment",   
    values_to = "Absorption")
l_SW1$Absorption <- as.numeric(l_SW1$Absorption)
l_SW1$Absorption <- ifelse(l_SW1$Absorption < 0, 0, l_SW1$Absorption)
l_SW1$lambda <- as.numeric(l_SW1$lambda)

l_SW2 <- l_SW2 %>%
  pivot_longer(
    cols = -1,             
    names_to = "Pigment",   
    values_to = "Absorption")
l_SW2$Absorption <- as.numeric(l_SW2$Absorption)
l_SW2$Absorption <- ifelse(l_SW2$Absorption < 0, 0, l_SW2$Absorption)
l_SW2$lambda <- as.numeric(l_SW2$lambda)

l_SW <- rbind(l_SW1, l_SW2) # (m2.mg-1)

# Plot 
# all_spectrums <- 
  ggplot(l_SW, aes(x = lambda, y = Absorption, color = Pigment)) +
  geom_line(size = 1) +
  theme_linedraw() +
  labs(
    x = colnames(df)[1],
    y = "Coefficient d'absorption (m2.mg-1)",
    title = "Pigment abs spectrums (Clementson  and Wojtasiewicz, 2019)")+
  ylim(0, 0.09)

# ggplotly(all_spectrums)

# Pigment selection
sub_pig <- c(
  "Antheraxanthin", 
  "B,B-carotene", 
  "B,e-carotene", 
  "chlorophyll c2", 
  "Chlorophyll-a", 
  "Diadinoxanthin", 
  "Diatoxanthin", 
  "Fucoxanthin", 
  "Lutein", 
  "Neoxanthin", 
  "Phaeophorbide-a", 
  "Phaeophytin-a",
  "Violaxanthin", 
  "Zeaxanthin")

# Retrieve absorption spectrums
abs_spectrum <- l_SW %>%
  filter(Pigment %in% sub_pig) 

# Rename
abs_spectrum <- abs_spectrum %>%
  mutate(
    Pigment = dplyr::recode(
      Pigment,
      `Phaeophytin-a` = "Pheophytin-a",
      `Phaeophorbide-a` = "Pheophorbide-a",
      `B,B-carotene` = "ββ",
      `B,e-carotene` = "βε",
      `chlorophyll c2` = "Chlc2",
      `Diadinoxanthin` = "Dd",
      `Diatoxanthin` = "Dt",
      `Lutein` = "L",
      `Neoxanthin` = "N",
      `Zeaxanthin` = "Z",
      `Antheraxanthin` = "A",
      `Violaxanthin` = "V"))

# Light spectrum range
abs_spectrum <- abs_spectrum %>% filter(lambda >= 400 & lambda <= 700) 

# Pigments used for the spectrum reconstruction
p_inter <- intersect(unique(abs_spectrum$Pigment), colnames(df_Qphar))

# Theorical spectral absorption a_mpb (dimensionless) ####

reconstruct_spectrum <- function(sample_row, abs_spectrum, p_inter) {
  
  abs_sample <- abs_spectrum %>%
    filter(Pigment %in% p_inter) %>% 
    mutate(Concentration = as.numeric(sample_row[ , Pigment]))
  
  a_rec <- abs_sample %>%
    mutate(AbsWeighted = Absorption * Concentration) %>%      # a_mpb calculation
    
    group_by(lambda) %>% 
    summarise(AbsorptionTotal = sum(AbsWeighted, na.rm = TRUE)) %>% 
    ungroup() %>%
    mutate(Sample = sample_row$File)
  
  return(a_rec)}

# Calculations of all spectrums
all_spectra_S <- map_df(1:nrow(df_Qphar_Surf), function(i) { 
  reconstruct_spectrum(df_Qphar_Surf[i, ], abs_spectrum, p_inter)}) 

# Plot of reconstitued uncorrected package-effect spectrums
all_spectra_S <- all_spectra_S %>%mutate(Condition = str_sub(Sample, 1, -3)) 

# a_mpb_un <- 
ggplot(all_spectra_S, aes(x = lambda, y = AbsorptionTotal, color = Condition, fill = Condition)) +

  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +

  stat_summary(
    fun.data = mean_sdl,       
    fun.args = list(mult = 1),  
    geom = "ribbon",
    alpha = 0.4,
    color = NA) +
  
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Surfacic absorption spectrum rec,un (dimessionless)",
    title = "Average spectrums without cell packing correction")

# ggplotly(a_mpb_un)

# Average DA WSB spectrum (dimensionless)
DA_WSB <- all_spectra_S[all_spectra_S$Condition == "DA_WSB", ]

ggplot(DA_WSB, aes(x = lambda, y = AbsorptionTotal, color = Condition, fill = Condition)) +

  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +

  stat_summary(
    fun.data = mean_sdl,       
    fun.args = list(mult = 1),  
    geom = "ribbon",
    alpha = 0.2,
    color = NA) +
  
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Surfacic absorption spectrum rec,un (dimessionless)",
    title = "Average WSB DA spectrums without cell packing correction")


# Incident LED panel spectrum ####

file_path_spectrum <- "1_Data_formatting/data/Qphar/1500PAR.xlsx" #spectrum

# Retrieve values
spectre <- read_excel(file_path_spectrum)
spectre <- spectre[-(1:54), ]

colnames(spectre)[colnames(spectre) == "Name"] <- "Wavelength"
spectre$Wavelength <- as.numeric(spectre$Wavelength)
colnames(spectre)[colnames(spectre) == "Value"] <- "Irradiance" # (W.m2)
spectre$Irradiance <- as.numeric(spectre$Irradiance)
spectre <- spectre %>% filter(`Wavelength`  >= 400 & `Wavelength`  <= 700) 

# W.m-2 -> µmol photons m-2 s-1 conversion (https://www.berthold.com/fr-fr/bio-analyses/base-de-connaissances/faq/comment-convertir-lirradiation-en-flux-de-photons/)
planck <- 6.62607015e-34    # J.s
celerite <- 2.99792458e8      # m.s
avogadro <- 6.02214076e23  # 1 mol photon = 6.02214076e23 photons

# PAR calculation at each waveslength
spectre <- spectre %>% mutate(Q_umol = ((`Irradiance` * (`Wavelength` * 1e-9)) / (planck * celerite) )    / (avogadro * 1e-6))

unique_lambda <- unique(spectre$Wavelength)

# Check the correspondence with the raw PPFD file value
raw_par <- sum(spectre$Q_umol, na.rm = TRUE) 
raw_par

# Incident spectrum plot 
total_PPFD <- sum(spectre$Q_umol, na.rm = TRUE) 

ggplot(spectre, aes(x = Wavelength, y = Q_umol)) +
  geom_line(size = 1) +
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Irradiance (µmol photons.m-2.s-1)",
    title = "Incident HL irradiance spectrum"
  ) +
  annotate(
    "text", 
    x = min(spectre$Wavelength),  
    y = max(spectre$Q_umol), 
    label = paste0("Total PPFD = ", round(total_PPFD, 1), " µmol photons m-2 s-1"),
    hjust = 0, vjust = 1,  
    size = 5)

# Cellular packing correction calculation (utilisation de df_Qphar_Vol) ####
#10.1029/2022JC018494 (Soja-wozniak et al., 2022)

# Eq. 1 (m-1)
all_spectra_V <- map_df(1:nrow(df_Qphar_Vol), function(i) { 
  reconstruct_spectrum(df_Qphar_Vol[i, ], abs_spectrum, p_inter)}) 

# Metadata formatting 
all_spectra_V <- all_spectra_V %>%mutate(Condition = str_sub(Sample, 1, -3)) 

# Average uncorrected DA WSB spectrum from volumic biofilm pigment concentration
a_rec_un <- all_spectra_V %>%filter(grepl("DA_WSB", Condition))

a_rec_un <- a_rec_un %>%
  group_by(lambda) %>%     
  summarise(abs_rec_un = mean(AbsorptionTotal, na.rm = TRUE)) %>%
  ungroup()

# Plot
ggplot(a_rec_un, aes(x = lambda, y = abs_rec_un)) +

  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +

  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "mean surfacic absorption spectrum rec,un (m-1)",
    title = "Average DA spectrums without cell packing correction")


# Eq. 6 Estimated MPB absorption at 675nm

# MPB absorbance at 675 nm (in vivo chla abs) (m-1)
a_rec_675 <- a_rec_un %>%
  summarise(a675 = approx(lambda, abs_rec_un, xout = 675)$y) %>% 
  pull(a675)

# Retrieve average chla concentration of DA WSB replicates (mg.m-3)
Ca_mgm3 <- mean(df_Qphar_Vol$`Chlorophyll-a`[df_Qphar_Vol$Treatment == "DA_WSB"])



# Average diatoms (P. strigosum) biovolume calculation

d <- "1_Data_formatting/data/Qphar/Pstrigosum_size.xlsx"
d <- read_excel(d)
d <- as.data.frame(t(colMeans(d, na.rm = TRUE)))

L <- d$Length # n = 30 (my measurments)
W <- d$Width # n = 30 (my measurments)
H1 <- L/11.18
H2 <- W/1.79 
H <- (H1+H2)/2 # Assuming a ratio (de tommasi et al., doi.org/10.1038/s41598-024-56206-y : L ave = (115+254)/2 ; W ave = (25+34)/2 ; H (13+20)/2 ; ratio L:H = 11.18 et W:H 1.79. Empirical data

# Diatom prism parallelogram-base calculation from doi.org/10.1046/j.1529-8817.1999.3520403.x (µm3)
V_diato <- 0.5*L*W*H 

# Spherical ESD biovolume conversion
r_eq <- (3 * V_diato / (4*pi))^(1/3) 
S_eq <- 4 * pi * r_eq^2 
S_V_diato <- S_eq / V_diato 

# Re-adjusting value for my large cells, extrapolation of microplankton Utiz et al., 2006 data

r_micro <- 35/2 #35 = ESD diameter. The specific absorption coefficient of Utiz et al. (2006) (0.009) is defined for microplankton (>20 µm ESD) and is therefore not directly applicable to my large diatoms. Because this value is associated with a broad size class (>20 µm) without further specification, it is necessary to approximate an appropriate ESD while remaining within model constraints.Assuming 0.009 corresponds to a given microplankton S/V ratio, I adjusted this coefficient proportionally to the S/V of my cells. Numerical tests showed that the ρ_root solution remains physically realistic only for cell sizes between ~29 and 40 µm; beyond 40 µm, ρ_root falls below 0.2, which would unrealistically underestimate package effects. For this reason, an intermediate reference size of 35 µm was selected, as it ensures convergence of the model while maintaining plausible photophysical properties.

V_micro <- 4/3 * pi * r_micro^3 
S_micro <- 4 * pi * r_micro^2
S_V_micro <- S_micro / V_micro

# Calculation of my coefficient (The value is less than 0.009, which makes sense because my cells were larger than typical microplankton)
a_star_diato <- 0.009 * (S_V_diato / S_V_micro) 

# Theorical absorption estimation of DA WSB at 675nm from real pigment concentrations and a_star_diato coefficient
a_estim_ph675 <- a_star_diato*Ca_mgm3*1 # (mg.m-3) multiplied by 1, assuming that my biovolume community was 100% dominated by large diatoms

# Eq. 7  Q*_a estimation at 675 nm
Qstar_est_675 <- a_estim_ph675 / a_rec_675 # (mg.m-2) 

# Eq. 8 at 675 nm
f_rho <- function(rho, Qstar) {
  (2/3) * Qstar * rho - 1 + 2 * (1 - (1 + rho) * exp(-rho)) / (rho^2)
}

lower <- 0.02 # upper and lower realistic range value
upper <- 1
fr_low <- f_rho(lower, Qstar_est_675)
fr_high <- f_rho(upper, Qstar_est_675)

if(fr_low * fr_high > 0) {
  warning("Error")
  rho_root <- NA_real_
} else {
  rho_root <- uniroot(function(x) f_rho(x, Qstar_est_675), c(lower, upper))$root
}


# Retrieve Chla absorption from clementson paper
a_star_chla_675 <- l_SW %>%
  filter(Pigment == "Chlorophyll-a") %>%
  summarise(a = approx(lambda, Absorption, xout = 675)$y) %>%
  pull(a)

# Eq. 9  Ci.d calculation at 675 nm
Ci_d <- rho_root / a_star_chla_675 # (mg.m-2)


# Eq. 5 Ci_d at 675 nm application at my whole reconstructed spectrum
a_rec_un <- a_rec_un %>%
  mutate(
    a_rec_star = abs_rec_un / Ca_mgm3,  #  (m²·mg⁻¹) Chla spectrum normalisation
    ρ = a_rec_star * Ci_d  )    # (unitless) Spectrum correction, where Ci_d respresent the effective absorptivity surfacic density of Chla unit at 675 nm


# Eq. 3 and 2 calculating Q*a at the whole spectrum
a_rec_un <- a_rec_un %>%
  mutate(
    Qa = 1 + (2 * exp(-`ρ`) / `ρ`) + 2 * ((exp(-`ρ`) - 1) / `ρ`^2), 
    Qastar = (3/2) * (Qa / ρ)) 


# Eq. 10 Package-effect correction

FL <- 1 # assuming that 100% of my MPB biovolume was concerned by the package effect

a_rec_un <- a_rec_un %>%
  mutate(abs_rec_corrected_pack = FL * abs_rec_un * Qastar + (1 - FL) * abs_rec_un) # (m-1) 

# Plot
spectre1 <- ggplot(a_rec_un, aes(x = lambda)) +
  geom_line(aes(y = abs_rec_un, color = "Original")) +
  geom_line(aes(y = abs_rec_corrected_pack, color = "Package-corrected")) +  
  labs(x = "Wavelength (nm)", y = "Reconstitued spectrum (m-1)", color = "Reconstituted spectrum", title = "") +
  theme_test()+
  theme(legend.position = "bottom")


# % loss calculation index at each wavelength
a_rec_un <- a_rec_un %>%
  mutate(decrease_package_effect = 100 - ((abs_rec_corrected_pack * 100) / abs_rec_un))


# Plot
spectre2 <- ggplot(a_rec_un, aes(x = lambda)) +
  geom_line(aes(y = decrease_package_effect)) +
  labs(x = "Wavelength (nm)", y = "Decrease pEffect (%)", title = "DA WSB") +
  theme_test()

spectre1+spectre2

# Retrieve % index loss 
spectrum_pCor_corrected <- data.frame(
  lambda = a_rec_un$lambda,                          
  decrease_package_effect = a_rec_un$decrease_package_effect /100) # (dimensionless)


# Values interpolation 
interp_cValues <- approx(
  x = spectrum_pCor_corrected$lambda,          
  y = spectrum_pCor_corrected$decrease_package_effect, 
  xout = unique_lambda,
  rule = 2
)$y

spectrum_pCor_interpolated <- data.frame(
  lambda = unique_lambda,
  decrease_package_effect = interp_cValues)

# Cellular packing correction application ####

# Retrieve uncorrected absorption spectrum
a_mpb <- all_spectra_S

a_mpb_interp <- a_mpb %>%
  group_by(Sample) %>%
  group_modify(~ {
    df <- .x
    
    if(nrow(df) == 0) {
      return(data.frame(lambda = unique_lambda,
                        AbsorptionTotal = rep(NA_real_, length(unique_lambda))))
    }
    
    agg <- df %>%
      group_by(lambda) %>%
      summarise(AbsorptionTotal = mean(AbsorptionTotal, na.rm = TRUE), .groups = "drop")
    
    if(nrow(agg) < 2) {
      return(data.frame(lambda = unique_lambda,
                        AbsorptionTotal = rep(NA_real_, length(unique_lambda))))
    }
    
    interp_vals <- approx(x = agg$lambda,
                          y = agg$AbsorptionTotal,
                          xout = unique_lambda,
                          rule = 2)$y
    
    data.frame(lambda = unique_lambda,
               AbsorptionTotal = interp_vals)
  }) %>%
  ungroup() %>%
  arrange(Sample, lambda)


# Metadata formatting
a_mpb_interp <- a_mpb_interp %>%mutate(Condition = str_sub(Sample, 1, -3)) 

# Spectrum subset
a_mpb_interp <- a_mpb_interp[grepl("DA_WSB", a_mpb_interp$Condition), ] 

# Plot
before_cor <- ggplot(a_mpb_interp, aes(x = lambda, y = AbsorptionTotal, color = Condition, fill = Condition)) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +
  
  stat_summary(
    fun.data = mean_sdl,       
    fun.args = list(mult = 1),  
    geom = "ribbon",
    alpha = 0.4,
    color = NA) +
  
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Surfacic absorption spectrum rec,un (dimensionless)",
    title = "Reconstructed spectrums without cell packing correction")



# Package effect loss correction at each wavelength
idx <- match(a_mpb_interp$lambda, spectrum_pCor_interpolated$lambda)

a_mpb_interp$AbsorptionTotal_corrected <- a_mpb_interp$AbsorptionTotal - a_mpb_interp$AbsorptionTotal * spectrum_pCor_interpolated$decrease_package_effect[idx]

# Plot
after_cor <- ggplot(a_mpb_interp, aes(x = lambda, y = AbsorptionTotal_corrected, color = Condition, fill = Condition)) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +
  
  stat_summary(
    fun.data = mean_sdl,       
    fun.args = list(mult = 1),  
    geom = "ribbon",
    alpha = 0.4,
    color = NA) +
  
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Surfacic absorption spectrum rec,pac (dimensionless)",
    title = "Reconstructed spectrums with cell packing correction")

before_cor+after_cor


# Qphar calculation ####

Qphar <- a_mpb_interp %>%
  inner_join(
    spectre,
    by = c("lambda" = "Wavelength"))


# Uncorrected Qphar calculation
Qphar <- Qphar %>%
  mutate(Qphar = Q_umol- (Q_umol * exp(-AbsorptionTotal))) 

# Pakage-effect corrected Qphar calculation
Qphar <- Qphar %>%
  mutate(Qphar_cor = Q_umol - (Q_umol * exp(-AbsorptionTotal_corrected)))

# Sum calculation for each sample
Qphar_sum <- Qphar %>%
  group_by(Sample) %>%
  summarise(Sum_Qphar = sum(Qphar, na.rm = TRUE))
Qphar_sum

Qphar_sum_cor <- Qphar %>%
  group_by(Sample) %>%
  summarise(Sum_Qphar_cor = sum(Qphar_cor, na.rm = TRUE))
Qphar_sum_cor

# Photonic energy calculation from Qphar_cor E= ∑(Photon flux λ) × hc/λ 
Ep <- (planck * celerite) / (unique_lambda* 1e-9) # J/photon

Ephoton <- data.frame( 
  lambda = unique_lambda,
  Ep = Ep)

Qphar$Eabs <- (Qphar$Qphar_cor * 1800 * (pi * (0.055/2)^2)) # 1800 = Time exposure (s) multiplied by Petri Dish (m2)

Qphar$Eabs <- (Qphar$Eabs * 1e-6 * avogadro) # µmol -> mol photons conversion mulitplied by avogadro to retrive tot absorbed photon by diatom suspension during the experimentation

# Df formatting
Qphar_photons <- Qphar

# Energy calculation at each wavelength from specific photon Energy
Qphar <- merge(Qphar, Ephoton, by = "lambda", all.x = TRUE)
Qphar$Eabs <- (Qphar$Eabs * Qphar$Ep) 

# Absorbed energy sum (J)
E_abs <- Qphar %>%
  group_by(Sample) %>%
  summarise(across(9, ~sum(.x, na.rm = TRUE)))

# Sum of the total photon number absorbed
n_Photons <- Qphar_photons %>%
  group_by(Sample) %>%
  summarise(across(9, ~sum(.x, na.rm = TRUE)))

n_Photons <- n_Photons %>% 
  rename(Pabs = Eabs) 

# Metadata formatting
E_abs <- E_abs %>%mutate(Condition = str_sub(Sample, 1, -3))

E_abs$Condition <- E_abs$Condition %>%
  gsub("DA_WSB", "Dark-Adapted Samples Without Sediment", .) 

# Table
light_table <- Qphar_sum %>%
  left_join(Qphar_sum_cor, by = "Sample")

# Merge data
light_table <- light_table %>%
  left_join(n_Photons, by = "Sample") %>%
  left_join(E_abs, by = "Sample")

# Mean sd calculations 
light_table <- light_table %>%
  group_by(Condition) %>%
  summarise(across(
    2:5,  
    list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd   = ~ sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(across(ends_with("_mean"), round, 2),
         across(ends_with("_sd"), round, 2)) %>%
  rowwise() %>%
  
  mutate(
    
    Pabs_mean = paste0(
      formatC(Pabs_mean, format = "e", digits = 2),   
      " (±",
      formatC(Pabs_sd, format = "e", digits = 2),
      ")"
    ),
    
    Sum_Qphar_mean = paste0(Sum_Qphar_mean, " (±", Sum_Qphar_sd, ")"),
    Sum_Qphar_cor_mean = paste0(Sum_Qphar_cor_mean, " (±", Sum_Qphar_cor_sd, ")"),
    Eabs_mean = paste0(Eabs_mean, " (±", Eabs_sd, ")")
  ) %>%
  
  select(Condition, ends_with("_mean")) %>%
  rename_with(~ gsub("_mean", "", .x), ends_with("_mean"))

# Add incident PPFD data
light_table <- light_table %>%
  mutate(
    PPFD = case_when(
      Condition == "Dark-Adapted Samples Without Sediment" ~ 1500,
      TRUE ~ NA_real_ )) %>%
  relocate(PPFD, .after = 1)

light_table <- light_table %>% rename(Qphar = Sum_Qphar)
light_table <- light_table %>% rename(`Qphar corrected` = Sum_Qphar_cor)

light_table %>%
  gt()%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left", weight = "bold"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgrey"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))