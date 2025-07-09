#Library packages ####
#install.packages("reshape2")
library(reshape2)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("cowplot")
library(cowplot)
#install.packages("patchwork")
library(patchwork)
#install.packages("ggrepel")
library(ggrepel)
#install.packages("readxl")
library(readxl)
#install.packages("sf")
library(sf)
#install.packages("stringr")
library(stringr)
#install.packages("gt")
library(gt)
#install.packages("tidyr")
library(tidyr)
#install.packages("vegan")
library(vegan)
#install.packages("dplyr")
library(dplyr)
#install.packages("rnaturalearth")
library(rnaturalearth)
#install.packages("here")
library(here)
#install.packages("webshot2")
library(webshot2)

# FIGURE 1 MPB DIV ####

  # MAP (A)

lat <- 47.875
lon <- -3.918

crs <- paste0("+proj=laea +lat_0=", lat, " +lon_0=", lon, " +datum=WGS84 +units=m +no_defs")

ctrys50m <- ne_countries(scale = 50, type = "countries", returnclass = "sf") %>%
  select(iso_a3, iso_n3, admin) %>%
  st_transform(crs = crs)

sphere <- st_graticule(ndiscr = 10000, margin = 10e-6) %>%
  st_transform(crs = crs) %>%
  st_convex_hull() %>%
  summarise(geometry = st_union(geometry))

concarneau <- data.frame(lon = lon, lat = lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = crs)

coords <- st_coordinates(concarneau)

marge <- 1e6 # Map limits 

xlim <- c(coords[1] - 7e5, coords[1] + 8.5e5)
ylim <- c(coords[2] - marge, coords[2] + marge)

ggplot() +
  geom_sf(data = sphere, fill = "#e2ebf7", alpha = 0.7) +
  geom_sf(data = ctrys50m, fill = "#ffffff", color = "black") + 
  geom_sf(data = concarneau, color = "red", size = 3) +
  coord_sf(crs = crs, xlim = xlim, ylim = ylim) +
  theme_bw()+
  theme(
    panel.grid.major = element_line(color = "black", size = 0.3)
  )+
  theme(
    axis.text.x = element_text(color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5))+
  plot_annotation(
    title = "(A)",
    caption = "") & theme(plot.caption = element_text(hjust = 0.5, size = 16, face = "bold"),
                                    plot.title = element_text(face="bold", size="20"))

  # BARPLOT (D)

df_18 <- psmelt(final_biofilm18S)

# Phylum group
df_agg_18 <- df_18 %>%
  group_by(Sample, Conditions, Biofilm_type, Class) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# Abundance means
df_agg_18 <- df_agg_18 %>%
  group_by(Conditions, Class) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = 'drop')

# Relative abundances
df_agg_18 <- df_agg_18 %>%
  group_by(Class, Conditions) %>%
  summarise(Mean_Abundance = sum(Mean_Abundance), .groups = 'drop') %>%
  group_by(Conditions) %>%
  mutate(Relative_Abundance = Mean_Abundance / sum(Mean_Abundance)*100)

# Threshold
df_agg_18 <- df_agg_18 %>%
  mutate(Class = ifelse(Relative_Abundance < 2, "Others", Class)) 

# Cleaning names
df_agg_18$Class <- gsub("\\s+", "", df_agg_18$Class)
df_agg_18 <- df_agg_18 %>%rename(Phylum = Class)

# df subset
df_WSB_18 <- subset(df_agg_18, grepl("^WSB_", Conditions))
df_SB_18 <- subset(df_agg_18, grepl("^SB_", Conditions))

# Colors settings
all_Phylum_18 <- union(unique(df_SB_18$Phylum), unique(df_WSB_18$Phylum))
color_Phylum_18 <- colorRampPalette(c("#0b090a","#1a759f", "#cae9ff", "#0a9396", "#a1cca5",  "#e9ecef", "#5fa8d3", "#708d81", "#005f73", "#83c5be", "grey")) (11)
Phylum_colors_18 <- data.frame(Phylum = all_Phylum_18, Color = color_Phylum_18)

# Phylum names
labels_italics <- setNames(
  lapply(all_Phylum_18, function(x) bquote(italic(.(x)))),
  all_Phylum_18)

# Plot SB
plot_SB_phylum_18 <- ggplot(df_SB_18, aes(x = Conditions, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.7) +
  theme_test() +
  labs(title = "SB", x = "", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = setNames(Phylum_colors_18$Color, Phylum_colors_18$Phylum),
                    labels = labels_italics)+
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom")+
  annotate("text", 
           x = 1, y = -4, 
           label = "VLL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 2, y = -4, 
           label = "HL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 3, y = -4, 
           label = "HP", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 4, y = -4, 
           label = "PL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")

# PLot WSB
plot_WSB_phylum_18 <- ggplot(df_WSB_18, aes(x = Conditions, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.7) +
  theme_test() +
  labs(title = "WSB", x = "", y = "") +
  theme(plot.title = element_text(face = "bold", size = 24),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = setNames(Phylum_colors_18$Color, Phylum_colors_18$Phylum))+
  annotate("text", 
           x = 1, y = -4, 
           label = "VLL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 2, y = -4, 
           label = "HL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 3, y = -4, 
           label = "HP", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 4, y = -4, 
           label = "PL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")

(plot_SB_phylum_18 + plot_WSB_phylum_18)+
  plot_annotation(
    title = "(D)",
    caption = "Conditions") & theme(plot.caption = element_text(hjust = 0.5, size = 16, face = "bold"),
                                    plot.title = element_text(face="bold", size="20"))

  # Pie chart (E)

df_WSB_18_pc <- df_WSB_18 %>%
  group_by(Phylum) %>%
  summarise(sum = sum(Mean_Abundance), .groups = 'drop') %>% 
  mutate(mean = 100 * sum / sum(sum))

# Label position 
df_WSB_18_pc <- df_WSB_18_pc %>%
  arrange(desc(Phylum)) %>%
  mutate(pos = cumsum(mean) - mean / 2) 

# Plot 1
ggplot(df_WSB_18_pc, aes(x = "", y = mean, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.9, position = "stack") +
  coord_polar(theta = "y", direction = -1) +  
  geom_label_repel(aes(label = paste0(round(mean, 2), " %"), y = pos),
                   color = "white", size = 10, fontface = "bold", 
                   nudge_x = 0.73,              
                   force = 0.6,                 
                   box.padding = 0.2,
                   show.legend = FALSE,
                   segment.color = "black") +
  theme_void() +
  labs(title = "(E)") +
  theme(plot.title = element_text(face = "bold", size = 16),
                                  legend.position = "none") +
  scale_fill_manual(values = setNames(Phylum_colors_18$Color, Phylum_colors_18$Phylum))

  # Pleurosigma domination, Pie chart n°2

dfeuk <-df_18

# Cleanning df
dfeuk[, 10:18] <- lapply(dfeuk[, 10:18], function(x) gsub("\\s+", "", x))

# Bacillariophyceae
dfeuk <- dfeuk %>% filter(Class %in% "Bacillariophyceae")

# Unidentified Bacillariophyceae correction (blastn ncbi)
correction_file <- here("data", "Metabarcoding_data", "unidentified_bacill_correction_blastn.xlsx")
correction_file <- read_excel(correction_file)
correction_file <- correction_file[, 1:2] #recup col 1 et 2
correction_file <- na.omit(correction_file) #suppr des taxons non identifiés par blastn

# df format
correction_file <- correction_file %>%
  separate(
    blastn_id_order, 
    into = c("Order", "Family", "Genus", "Species"),
    sep = ";",
    remove = TRUE)

# merge & correction
dfeuk <- dfeuk %>%
  left_join(correction_file, by = c("OTU" = "ID_ASV"))

dfeuk <- dfeuk %>%
  mutate(
    Order = coalesce(Order.y, Order.x),  
    Family = coalesce(Family.y, Family.x),
    Genus = coalesce(Genus.y, Genus.x),
    Species = coalesce(Species.y, Species.x)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Subset
dfeuk <- subset(dfeuk, grepl("^WSB_", Conditions))

# Species
dfeuk_sp <- dfeuk[, c("Species", "Abundance")]
dfeuk_sp$Species <- sub("_n.*", "", dfeuk_sp$Species) # Tag suppr

# Calcul by species and % transformation
dfeuk_sp <- dfeuk_sp %>%
  group_by(Species) %>%
  summarise(sum = sum(Abundance), .groups = 'drop') %>%
  mutate(mean = 100 * sum / sum(sum))

# Name format
dfeuk_sp <- dfeuk_sp %>% mutate(Species = gsub("_", " ", Species))

# Threshold 
dfeuk_sp <- dfeuk_sp %>%mutate(Species = ifelse(mean < 2, "Others", Species))

# merge other
dfeuk_sp <- dfeuk_sp %>%
  group_by(Species) %>%
  summarise(sum = sum(mean), .groups = 'drop') %>%
  mutate(mean = 100 * sum / sum(sum))

# Color settings
ord_sp_bacill <- union(unique(dfeuk_sp$Species), unique(dfeuk_sp$Species))
col_bacill <- c(
  "IncertaeSedis" = "grey", 
  "Other" = "#0b090a", 
  "Pleurosigma strigosum" = "#31572c")
ord_bacill <- data.frame(Species = ord_sp_bacill, Color = col_bacill) 

# Label position
dfeuk_sp <- dfeuk_sp %>%
  arrange(desc(Species)) %>%
  mutate(pos = cumsum(mean) - mean / 2)

# Species names
labels_italics <- setNames(
  lapply(ord_bacill$Species, function(x) bquote(italic(.(x)))),
  ord_bacill$Species)

ggplot(dfeuk_sp, aes(x = "", y = mean, fill = Species)) +
  geom_bar(stat = "identity", width = 1, position = "stack", color = "black", size = 1) +#0.7 legend et 1 borders
  coord_polar(theta = "y", start = pi) +
  geom_label_repel(aes(label = paste0(round(mean, 2), " %"), y = pos),
                   color = "white", size = 10, fontface = "bold",
                   nudge_x = 0.7,              
                   force = 0.5,                 
                   box.padding = 0.2,
                   show.legend = FALSE,
                   segment.color = "black") +
  theme_void() +
  labs(title = "") +
  theme(legend.position = "right") +
  scale_fill_manual(values = setNames(ord_bacill$Color, ord_bacill$Species),
                    labels = labels_italics)

# FIGURE 2 MIGRATION ####

y_limits <- c(2000, 18900)

# Plot HL
plot_hl <- ggplot(File_HL_data, aes(x = minutes, y = `F0`, group = Conditions, color = Conditions)) + 
  annotate("rect", xmin = 13, xmax = 43, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(group = Conditions), method = "gam", se = TRUE, formula = y ~ s(x), span = 0.4) +
  geom_ribbon(aes(ymin = ..ymin.., ymax = ..ymax.., group = Conditions, fill = Conditions),
              stat = "smooth", method = "gam", se = TRUE, formula = y ~ s(x), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#E7B75F", "#324851")) + 
  scale_fill_manual(values = c("#E7B75F", "#324851")) +  
  labs(x = "", y = "Surface microalgal biomass (F0 a.u.)", 
       title = "", 
       color = "Conditions", fill = "Conditions") + 
  theme_test()+
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),  
        axis.text.x = element_text(size = 18),                
        axis.text.y = element_text(size = 18)) +
  ylim(y_limits)

# Plot HP
plot_hp <- ggplot(File_HP_data, aes(x = minutes, y = `F0`, group = Conditions, color = Conditions)) + 
  annotate("rect", xmin = 25, xmax = 32.3, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(group = Conditions), method = "gam", se = TRUE, formula = y ~ s(x), span = 0.4) + 
  geom_ribbon(aes(ymin = ..ymin.., ymax = ..ymax.., group = Conditions, fill = Conditions),
              stat = "smooth", method = "gam", se = TRUE, formula = y ~ s(x), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#5C88C4", "#324851")) + 
  scale_fill_manual(values = c("#5C88C4", "#324851")) +  
  labs(x = "Time (Minutes)", y = "", 
       title = "", 
       color = "Conditions", fill = "Conditions") + 
  theme_test()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),  
        axis.text.x = element_text(size = 18),                
        axis.text.y = element_blank()) +
  theme(legend.position = "bottom") +
  ylim(y_limits)

# Plot PL
plot_pl <- ggplot(File_PL_data, aes(x = minutes, y = `F0`, group = Conditions, color = Conditions)) + 
  annotate("rect", xmin = 20, xmax = 31, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(group = Conditions), method = "gam", se = TRUE, formula = y ~ s(x), span = 0.4) + 
  geom_ribbon(aes(ymin = ..ymin.., ymax = ..ymax.., group = Conditions, fill = Conditions),
              stat = "smooth", method = "gam", se = TRUE, formula = y ~ s(x), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#937DC2", "#324851")) + 
  scale_fill_manual(values = c("#937DC2", "#324851")) +  
  labs(x = "", y = "", 
       title = "", 
       color = "Conditions", fill = "Conditions") + 
  theme_test()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),  
        axis.text.x = element_text(size = 18),                
        axis.text.y = element_blank()) +
  theme(legend.position = "bottom") +
  ylim(y_limits)

# Final plot
(plot_hl | plot_hp | plot_pl) 

# FIGURE 3 HEATMAP ####

df_heatmap <- df_Pparameters %>%select(Time, Name, Alpha, `rETRm obs`, "Fq'/Fm' VLL", NPQ1000, YNPQ1000, `YNOm obs`)

# Mobility subset
df_heatmap_wsb <- subset(df_heatmap, grepl("^WSB_", Name))
df_heatmap_sb <- subset(df_heatmap, grepl("^SB_", Name))

# Differences After-Before
df_heatmap_sb_diff <- df_heatmap_sb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = Alpha[Time == "After"] - Alpha[Time == "Before"],
    `Fq'/Fm' VLL` = `Fq'/Fm' VLL`[Time == "After"] - `Fq'/Fm' VLL`[Time == "Before"],
    `rETRm obs` = `rETRm obs`[Time == "After"] - `rETRm obs`[Time == "Before"],
    `NPQ 1000` = NPQ1000[Time == "After"] - NPQ1000[Time == "Before"],
    `Y(NPQ) 1000` = YNPQ1000[Time == "After"] - YNPQ1000[Time == "Before"],
    `Y(NO)m obs` = `YNOm obs`[Time == "After"] - `YNOm obs`[Time == "Before"])

df_heatmap_wsb_diff <- df_heatmap_wsb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = Alpha[Time == "After"] - Alpha[Time == "Before"],
    `Fq'/Fm' VLL` = `Fq'/Fm' VLL`[Time == "After"] - `Fq'/Fm' VLL`[Time == "Before"],
    `rETRm obs` = `rETRm obs`[Time == "After"] - `rETRm obs`[Time == "Before"],
    `NPQ 1000` = NPQ1000[Time == "After"] - NPQ1000[Time == "Before"],
    `Y(NPQ) 1000` = YNPQ1000[Time == "After"] - YNPQ1000[Time == "Before"],
    `Y(NO)m obs` = `YNOm obs`[Time == "After"] - `YNOm obs`[Time == "Before"])

# SB comparaison groups by day
SB_VLL_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_VLL_1", "SB_VLL_2", "SB_VLL_3")) # J1 experiment
SB_HL_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_HL_1", "SB_HL_2", "SB_HL_3")) # J1 experiment
SB_VLL_4_6 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_VLL_4", "SB_VLL_5", "SB_VLL_6")) # J2 experiment
SB_PL_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_PL_1", "SB_PL_2", "SB_PL_3")) # J2 experiment
SB_VLL_7_9 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_VLL_7", "SB_VLL_8", "SB_VLL_9")) # J3 experiment
SB_HP_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_HP_1", "SB_HP_2", "SB_HP_3")) # J3 experiment

# WSB comparaison groups by day
WSB_VLL_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_VLL_1", "WSB_VLL_2", "WSB_VLL_3"))
WSB_HL_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_HL_1", "WSB_HL_2", "WSB_HL_3"))
WSB_VLL_4_6 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_VLL_4", "WSB_VLL_5", "WSB_VLL_6"))
WSB_PL_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_PL_1", "WSB_PL_2", "WSB_PL_3"))
WSB_VLL_7_9 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_VLL_7", "WSB_VLL_8", "WSB_VLL_9"))
WSB_HP_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_HP_1", "WSB_HP_2", "WSB_HP_3"))

# P-values loop function
perform_t_tests_pv_SB <- function(group1, group2) { # For SB group
  p_values <- c()
  for (param in colnames(df_heatmap_sb_diff)[2:7]) {
    t_test <- t.test(group1[[param]], group2[[param]])
    p_values <- c(p_values, t_test$p.value)
  }
  return(p_values)
}

perform_t_tests_pv_WSB <- function(group1, group2) { # For WSB group
  p_values <- c()
  for (param in colnames(df_heatmap_wsb_diff)[2:7]) {
    t_test <- t.test(group1[[param]], group2[[param]])
    p_values <- c(p_values, t_test$p.value)
  }
  return(p_values)
}

# P-values calculations stress vs VLL
p_values_HL_SB <- perform_t_tests_pv_SB(SB_VLL_1_3, SB_HL_1_3)
p_values_PL_SB <- perform_t_tests_pv_SB(SB_VLL_4_6, SB_PL_1_3)
p_values_HP_SB <- perform_t_tests_pv_SB(SB_VLL_7_9, SB_HP_1_3)

p_values_HL_WSB <- perform_t_tests_pv_WSB(WSB_VLL_1_3, WSB_HL_1_3)
p_values_PL_WSB <- perform_t_tests_pv_WSB(WSB_VLL_4_6, WSB_PL_1_3)
p_values_HP_WSB <- perform_t_tests_pv_WSB(WSB_VLL_7_9, WSB_HP_1_3)

# P-values df
df_heatmap_diff_pv <- data.frame(
  "HL" = p_values_HL_SB,
  "PL" = p_values_PL_SB,
  "HP" = p_values_HP_SB,
  "HL." = p_values_HL_WSB,
  "PL." = p_values_PL_WSB,
  "HP." = p_values_HP_WSB)

# Add photosynthetic parameters names metadata
df_heatmap_diff_pv <- t(df_heatmap_diff_pv)
colnames(df_heatmap_diff_pv) <- colnames(df_heatmap_sb_diff)[2:7]

# Df formatting
df_heatmap_diff_pv_long <- melt(as.matrix(df_heatmap_diff_pv))
colnames(df_heatmap_diff_pv_long) <- c("Condition", "Variable", "Value")

# Add p-values categories
df_heatmap_diff_pv_long$Categorie <- cut(df_heatmap_diff_pv_long$Value,breaks = c(-Inf, 0.05, 0.1, Inf),labels = c("< 0.05", "0.05 < > 0.1", "> 0.1"),right = FALSE)

# % differences After-Before
df_heatmap_sb_diff_perc <- df_heatmap_sb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = (((Alpha[Time == "After"]*100) / Alpha[Time == "Before"])-100),
    `Fq'/Fm' VLL` = (((`Fq'/Fm' VLL`[Time == "After"]*100) / `Fq'/Fm' VLL`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

df_heatmap_wsb_diff_perc <- df_heatmap_wsb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = (((Alpha[Time == "After"]*100) / Alpha[Time == "Before"])-100),
    `Fq'/Fm' VLL` = (((`Fq'/Fm' VLL`[Time == "After"]*100) / `Fq'/Fm' VLL`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

# Mean % loop function
perform_mean_diff_sb <- function(group1, group2) {
  mean_diff <- c()
  for (param in colnames(df_heatmap_sb_diff_perc)[2:7]) {
    diff <- mean(group2[[param]]) - mean(group1[[param]])
    mean_diff <- c(mean_diff, diff)
  }
  return(mean_diff)
}

perform_mean_diff_wsb <- function(group1, group2) {
  mean_diff <- c()
  for (param in colnames(df_heatmap_wsb_diff_perc)[2:7]) {
    diff <- mean(group2[[param]]) - mean(group1[[param]])
    mean_diff <- c(mean_diff, diff)
  }
  return(mean_diff)
}

# Groups
SB_VLL_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_VLL_1", "SB_VLL_2", "SB_VLL_3")) # J1 experiment
SB_HL_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_HL_1", "SB_HL_2", "SB_HL_3")) # J1 experiment
SB_VLL_4_6_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_VLL_4", "SB_VLL_5", "SB_VLL_6")) # J2 experiment
SB_PL_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_PL_1", "SB_PL_2", "SB_PL_3")) # J2 experiment
SB_VLL_7_9_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_VLL_7", "SB_VLL_8", "SB_VLL_9")) # J3 experiment
SB_HP_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_HP_1", "SB_HP_2", "SB_HP_3")) # J3 experiment

WSB_VLL_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_VLL_1", "WSB_VLL_2", "WSB_VLL_3"))
WSB_HL_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_HL_1", "WSB_HL_2", "WSB_HL_3"))
WSB_VLL_4_6_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_VLL_4", "WSB_VLL_5", "WSB_VLL_6"))
WSB_PL_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_PL_1", "WSB_PL_2", "WSB_PL_3"))
WSB_VLL_7_9_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_VLL_7", "WSB_VLL_8", "WSB_VLL_9"))
WSB_HP_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_HP_1", "WSB_HP_2", "WSB_HP_3"))

# Means calculations stress vs VLL
mean_diff_HL_SB_perc <- perform_mean_diff_sb(SB_VLL_1_3_perc, SB_HL_1_3_perc)
mean_diff_PL_SB_perc <- perform_mean_diff_sb(SB_VLL_4_6_perc, SB_PL_1_3_perc)
mean_diff_HP_SB_perc <- perform_mean_diff_sb(SB_VLL_7_9_perc, SB_HP_1_3_perc)

mean_diff_HL_WSB_perc <- perform_mean_diff_wsb(WSB_VLL_1_3_perc, WSB_HL_1_3_perc)
mean_diff_PL_WSB_perc <- perform_mean_diff_wsb(WSB_VLL_4_6_perc, WSB_PL_1_3_perc)
mean_diff_HP_WSB_perc <- perform_mean_diff_wsb(WSB_VLL_7_9_perc, WSB_HP_1_3_perc)

# Means df
df_heatmap_mean_diff_perc <- data.frame(
  "HL" = mean_diff_HL_SB_perc,
  "PL" = mean_diff_PL_SB_perc,
  "HP" = mean_diff_HP_SB_perc,
  "HL." = mean_diff_HL_WSB_perc,
  "PL." = mean_diff_PL_WSB_perc,
  "HP." = mean_diff_HP_WSB_perc)

# Add photosynthetic parameters names metadata
df_heatmap_mean_diff_perc <- t(df_heatmap_mean_diff_perc)
colnames(df_heatmap_mean_diff_perc) <- colnames(df_heatmap_sb_diff)[2:7]

# Df formatting
df_heatmap_mean_diff_perc_long <- melt(as.matrix(df_heatmap_mean_diff_perc))
colnames(df_heatmap_mean_diff_perc_long) <- c("Condition", "Variable", "Perc")

# Add arrows pictograms and colors
df_heatmap_mean_diff_perc_long$Arrow <- ifelse(df_heatmap_mean_diff_perc_long$Perc > 0, "↑", "↓")
df_heatmap_mean_diff_perc_long$Arrow_Color <- ifelse(df_heatmap_mean_diff_perc_long$Arrow == "↓", "firebrick", "black")

# Add 'Perc' % value column
df_heatmap_diff_pv_long <- cbind(df_heatmap_diff_pv_long, df_heatmap_mean_diff_perc_long[c("Arrow", "Arrow_Color", "Perc")])

# Converting 'Perc' (Absolute values)
df_heatmap_diff_pv_long$Perc <- abs(df_heatmap_diff_pv_long$Perc)

# Plot graph
ggplot(df_heatmap_diff_pv_long, aes(x = Condition, y = factor(Variable, levels = rev(unique(Variable))), fill = Categorie)) +
  geom_tile(color = "black", size = 0.7) + 
  scale_fill_manual(name = "P-values Welch Two Sample t-test", 
                    values = c("< 0.05" ="#6096ba", "0.05 < > 0.1" ="#a3cef1", "> 0.1"="white")) + 
  labs(title = "", x = "", y = "") +
  #geom_text(data = subset(df_heatmap_diff_pv_long, Categorie %in% c("0.05 < > 0.1", "< 0.05", "> 0.1")), 
            #aes(label = round(Perc, 2)), color = "black", size = 6,  hjust = 0.5, vjust = 0.3, fontface = "bold")+
  geom_text(data = subset(df_heatmap_diff_pv_long, Categorie != "> 0.1"), 
            aes(label = Arrow, color = Arrow_Color), size = 12,  hjust = 0.5, vjust = 0.4) +
  scale_color_identity()+
  theme_test()+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 16),
        axis.text.y = element_text(color = "black", face = "bold", size = 16),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  annotate("segment", 
           x = 0.51, xend = 3.48,  
           y = 6.8, yend = 6.8,  
           colour = "black", 
           size = 1) +
  annotate("segment", 
           x = 3.51, xend = 6.5,  
           y = 6.8, yend = 6.8,  
           colour = "black", 
           size = 1) +
  annotate("text", 
           x = 2, y = 7.3, 
           label = "SB", 
           size = 10, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 5, y = 7.3, 
           label = "WSB", 
           size = 10, 
           fontface = "bold", 
           colour = "black")+
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.37)))

# FIGURE 4 PPARAMETERS ####

# Retrieve parameters
df_long_Pp <- df_Pparameters %>%
  mutate(RowID = row_number()) %>%
  pivot_longer(cols = c(23:29, 37:50),
               names_to = "Pparameter",
               values_to = "Value") %>%
  select(1:2, Pparameter, Value)

# Df formatting
df_long_Pp <- df_long_Pp %>%
  mutate(Irradiance = case_when(
    grepl("_l1", Pparameter) ~ "10",
    grepl("_l2", Pparameter) ~ "20",
    grepl("_l3", Pparameter) ~ "50",
    grepl("_l4", Pparameter) ~ "100",
    grepl("_l5", Pparameter) ~ "300",
    grepl("_l6", Pparameter) ~ "500",
    grepl("_l7", Pparameter) ~ "1000"))

df_long_Pp <- df_long_Pp %>%
  mutate(Pparameter = case_when(
    grepl("^qy", Pparameter)   ~ "Fq'/Fm'",
    grepl("^YNO", Pparameter)  ~ "YNOm obs",
    grepl("^YNPQ", Pparameter) ~ "YNPQ1000"))

df_long_Pp <- df_long_Pp %>% filter(str_starts(Name, "WSB"))

df_long_Pp <- df_long_Pp %>% mutate(Stress = str_extract(Name, "(?<=_)[^_]+(?=_)"))

df_long_Pp <- df_long_Pp %>% mutate(Stress = paste(Stress, Time, sep = " "))

# Set Irradiance levels
df_long_Pp <- df_long_Pp %>%mutate(Irradiance = factor(Irradiance, levels = c(10, 20, 50, 100, 300, 500, 1000), ordered = TRUE))

# Set colors
Pp_colors <- c(
  "VLL Before" = "#CBC9C9",
  "VLL After" = "#040404",
  
  "HL Before" = "#EADD95",
  "HL After"  = "#E7B75F",
  
  "HP Before" = "#94B5FD",
  "HP After"  = "#5C88C4",
  
  "PL After"  = "#937DC2",
  "PL Before" = "#DAB7FE")

# Plot Fq'/Fm'
df_fqfm <- df_long_Pp %>% filter(Pparameter == "Fq'/Fm'")

aPp <- ggplot(df_fqfm, aes(x = Irradiance, y = Value, color = Stress, group = Stress)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  geom_vline(xintercept = "100", linetype = "dashed", color = "grey50")+
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "Fq'/Fm' (unitless)",
       color = "Stress") +
  theme_test() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# Plot YNPQ1000
df_ynpq <- df_long_Pp %>% filter(Pparameter == "YNPQ1000")

bPp <- ggplot(df_ynpq, aes(x = Irradiance, y = Value, color = Stress, group = Stress)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  geom_vline(xintercept = "100", linetype = "dashed", color = "grey50")+
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "Irradiance (µmol photons.m-2.s-1)",
       y = "Y(NPQ) (unitless)",
       color = "Stress") +
  theme_test() +
  theme(axis.title.x = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "bottom")

# Plot YNOm obs
df_yno <- df_long_Pp %>% filter(Pparameter == "YNOm obs")

cPp <- ggplot(df_yno, aes(x = Irradiance, y = Value, color = Stress, group = Stress)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  geom_vline(xintercept = "100", linetype = "dashed", color = "grey50")+
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "Y(NO) (unitless)",
       color = "Stress") +
  theme_test() +
  theme( legend.position = "none",
         axis.title.x = element_text(size = 18, face = "bold"), 
         axis.title.y = element_text(size = 18, face = "bold"),  
         axis.text.x = element_text(size = 14),
         axis.text.y = element_text(size = 14))

(aPp | bPp | cPp)
# FIGURE 5 PIGMENT BIPLOT ####

pigments <- df_pigment_perc %>% select(5:ncol(df_pigment_perc))
pca_pigments <- prcomp(pigments, scale = TRUE) # PCA

data <- cbind(df_pigment_perc, pca_pigments$x) # Add PCA results

n<-45 # Scale

source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R") # Convex function


ggplot(data, aes(x = -PC1, y = -PC2, color = Conditions, shape = Groups)) +
  geom_point(size = 3) +
  
  scale_color_manual(values = c("VLL" = "#324851", "HP" = "#5C88C4", "HL"="#E7B800", "PL"="#937DC2")) +
  scale_fill_manual(values = c("VLL" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2"))+
  geom_convexhull(aes(fill = Conditions), alpha=0.6)+
  stat_ellipse(aes(fill = Groups), type = "norm", level = 0.95, alpha = 0.2, geom = "polygon", linetype = "dashed", color = "black") + 
  annotate("segment",
           x = 0, y = 0,
           xend = -pca_pigments$rotation[,1] * n,
           yend = -pca_pigments$rotation[,2] * n,
           arrow = arrow(length = unit(0.2, "cm")),
           color = "darkgrey",
           linewidth = 0.2) +
  geom_text_repel(
    data = as.data.frame(pca_pigments$rotation),
    aes(x = -PC1 * n, y = -PC2 * n, label = colnames(pigments)),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = 100) +
  labs(
    x = paste0("PC 1 (", round(100 * pca_pigments$sdev[1]^2 / sum(pca_pigments$sdev^2), 1), "%)"),
    y = paste0("PC 2 (", round(100 * pca_pigments$sdev[2]^2 / sum(pca_pigments$sdev^2), 1), "%)")) +
  guides(
    color = guide_legend(override.aes = list(shape = NA)),
    shape = guide_legend(),
    fill = guide_legend(override.aes = list(shape = c(19, 17), color = "black", linetype = "dashed"))) +
  theme_test() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "bottom")


# Anosim, diff between SB and WSB in their pigment compositions
coord_axes <- pca_pigments$x[,1]
group_fac <- df_pigment_perc$Groups # metadata factor

coord_axes <- data.frame(coord_axes, group = group_fac)

# Euclidean distance matrix
distance_matrix_pigment <- dist(coord_axes$coord_axes, method = "euclidean")
anosim(distance_matrix_pigment, coord_axes$group, permutations = 9999)

# Simper
pigments_data <- df_pigment_perc[, -c(1:4)]   
grouping_var <- df_pigment_perc$Groups            

simper_result <- simper(pigments_data, grouping_var, permutations = 9999)

simper_df <- bind_rows(
  lapply(simper_result, function(x) {
    df <- as.data.frame(x)
    df$Pigments <- rownames(df)
    df[, c("Pigments", "average", "sd", "p")]
  }))

simper_df %>%
  gt() %>%
  tab_header(title = "SIMPER results – SB vs WSB contrast") %>%
  cols_label(
    Pigments = "Pigments",
    average = "Mean difference",
    sd = "SD",
    p = "p-value") %>%
  fmt_number(columns = c(average, sd, p), decimals = 4) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = list(
      cells_body(columns = everything()),
      cells_column_labels(columns = everything()),
      cells_title("title")))

# FIGURE 6 BARPLOT RATIOS ####
gdf_pigment_ratio <- gdf_pigment

# Ratio calculations
gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(rChla_deriv = (`Chla-derivatives` / (`Chla-derivatives` + Chla + Pheopigments)) * 100)

gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(rPheo = (Pheopigments / (`Chla-derivatives` + Chla + Pheopigments)) * 100)

gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(`D-DES` = (Dt / (Dd + Dt) * 100))

gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(`Z-DES` = (Z / (Z + A + V) * 100))

# Df barplot
barplot_ratios <- as.data.frame(gdf_pigment_ratio[, -c(5:34)]) 
barplot_ratios$Conditions <- factor(barplot_ratios$Conditions, levels = c("PL", "HP", "HL", "VLL"))

# Standard deviation and mean calculations
barplot_ratios_rChla_deriv <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(rChla_deriv),
    sd = sd(rChla_deriv),
    se = sd(rChla_deriv) / sqrt(n()))
barplot_ratios_rChla_deriv$Conditions <- factor(barplot_ratios_rChla_deriv$Conditions, levels = c("PL", "HP", "HL", "VLL"))

barplot_ratios_rPheo <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(rPheo),
    sd = sd(rPheo),
    se = sd(rPheo) / sqrt(n()))
barplot_ratios_rPheo$Conditions <- factor(barplot_ratios_rPheo$Conditions, levels = c("PL", "HP", "HL", "VLL"))

barplot_ratios_D_DES <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(`D-DES`),
    sd = sd(`D-DES`),
    se = sd(`D-DES`) / sqrt(n()))
barplot_ratios_D_DES$Conditions <- factor(barplot_ratios_D_DES$Conditions, levels = c("PL", "HP", "HL", "VLL"))

barplot_ratios_Z_DES <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(`Z-DES`),
    sd = sd(`Z-DES`),
    se = sd(`Z-DES`) / sqrt(n()))
barplot_ratios_Z_DES$Conditions <- factor(barplot_ratios_Z_DES$Conditions, levels = c("PL", "HP", "HL", "VLL"))

# Reverse values
barplot_ratios <- barplot_ratios %>%
  mutate(across(5:8, ~ ifelse(Groups == "SB", -.x, .x)))

barplot_ratios_rChla_deriv <- barplot_ratios_rChla_deriv %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

barplot_ratios_rPheo <- barplot_ratios_rPheo %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

barplot_ratios_D_DES <- barplot_ratios_D_DES %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

barplot_ratios_Z_DES <- barplot_ratios_Z_DES %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

# Colors
barplot_ratios_colors <- c("VLL" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2")

# Chla_deriv plot
a <- ggplot() +
  geom_bar(data = barplot_ratios_rChla_deriv, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.7) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = rChla_deriv, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_rChla_deriv,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = barplot_ratios_colors, breaks = c("VLL", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(face = "bold", size = 16),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none")+
  
  scale_y_continuous(breaks = seq(-30, 30, 10),
                     labels = abs(seq(-30, 30, 10))) +
  
  labs(title = "")+
  ylab("Mean Chla-derivatives normalized to Chla (%)")+
  theme(plot.title = element_text(face = "bold"))+
  
  annotate("text", 
           x =5, y = -20,
           label = "(A) Chla-derivatives", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x =5.5, y = -26,
           label = "", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x = 4, y = 25.2,
           label = "", 
           fontface = "bold", 
           size = 5)

# Pheo plot
b <- ggplot() +
  geom_bar(data = barplot_ratios_rPheo, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.7) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = rPheo, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_rPheo,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = barplot_ratios_colors, breaks = c("VLL", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(face = "bold", size = 16),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none")+
  
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = abs(seq(-10, 10, 5))) +
  
  labs(title = "")+
  ylab("Mean Pheopigments normalized to Chla (%)")+
  theme(plot.title = element_text(face = "bold"))+
  
  annotate("text", 
           x =5, y = -7,
           label = "(B) Pheopigments", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x =5.5, y = -9,
           label = "", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x = 4, y = 8.7,
           label = "", 
           fontface = "bold", 
           size = 5)
  

# D-DES plot
c <- ggplot() +
  geom_bar(data = barplot_ratios_D_DES, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.7) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = `D-DES`, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_D_DES,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = barplot_ratios_colors, breaks = c("VLL", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(face = "bold", size = 16),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none")+
  
  scale_y_continuous(breaks = seq(-50, 50, 10),
                     labels = abs(seq(-50, 50, 10))) +
  
  labs(title = "")+
  ylab("Mean D-DES (%)")+
  theme(plot.title = element_text(face = "bold"))+
  
  annotate("text", 
           x =5, y = -28,
           label = "(C) Diatoxanthin de-epoxidation state", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x = 5.5, y = -50,
           label = "", 
           fontface = "bold", 
           size = 5)
  
# Z-DES
d <- ggplot() +
  geom_bar(data = barplot_ratios_Z_DES, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.7) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = `Z-DES`, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_Z_DES,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = barplot_ratios_colors, breaks = c("VLL", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(face = "bold", size = 16),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom")+
  
  scale_y_continuous(breaks = seq(-40, 40, 10),
                     labels = abs(seq(-40, 40, 10))) +
  
  labs(title = "")+
  ylab("Mean Z-DES (%)")+
  theme(plot.title = element_text(face = "bold"))+
  
  annotate("text", 
           x =5, y = -20,
           label = "(D) Zeaxanthin de-epoxidation state", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x =5.5, y = 27,
           label = " ", 
           fontface = "bold", 
           size = 5)+
  annotate("text", 
           x = 4, y = 33.5,
           label = "", 
           fontface = "bold", 
           size = 5)

(a | c) /
(b | d)

# P-values calculations
wilcox.test(`Z-DES` ~ Conditions, data = barplot_ratios, 
            subset = Groups == "SB" & Conditions %in% c("VLL", "PL"),
            exact = TRUE)
# TABLE 1 ####

# HL values
hl_var <- File_HL_data %>%
  select(-c(1:4)) %>%
  filter((minutes >= 11 & minutes <= 13) | (minutes >= 44 & minutes <= 47)) %>%  
  mutate(Time = case_when(
    minutes >= 11 & minutes <= 13 ~ "Before",
    minutes >= 44 & minutes <= 47 ~ "After")) # Select before and after values

# F0 values mean before and after
hl_sum <- hl_var %>%
  group_by(Conditions, Time) %>%
  summarise(
    mean = mean(`F0`, na.rm = TRUE),
    sd = sd(`F0`, na.rm = TRUE),
    .groups = "drop")

# Format
hl_sum_wide <- hl_sum %>%
  pivot_wider(names_from = Time, values_from = c(mean, sd))

hl_sum <- hl_sum_wide %>%
  mutate(
    diff_After_Before = mean_After - mean_Before,
    percent_change = (diff_After_Before / mean_Before) * 100,
    diff_After_Before = sprintf("%+.0f", diff_After_Before),
    percent_change = sprintf("%+.0f%%", percent_change),
    Before = sprintf("%.0f (%.0f)", mean_Before, sd_Before),
    After = sprintf("%.0f (%.0f)", mean_After, sd_After)
  ) %>%
  select(Conditions, Before, After, diff_After_Before, percent_change)

# Add VLL (HL) tag
hl_sum <- hl_sum %>%
  mutate(Conditions = if_else(str_detect(Conditions, "VLL"), 
                          paste0(Conditions, " (HL)"), 
                          Conditions))

# HP values
hp_var <- File_HP_data %>%
  select(-c(1:4)) %>%
  filter(minutes >= 24 & minutes <= 33) %>%
  mutate(Time = case_when(
    minutes >= 24 & minutes <= 26 ~ "Before",
    minutes >= 30 & minutes <= 33 ~ "After"))

hp_sum <- hp_var %>%
  group_by(Conditions, Time) %>%
  summarise(
    mean = mean(`F0`, na.rm = TRUE),
    sd = sd(`F0`, na.rm = TRUE),
    .groups = "drop")

hp_sum_wide <- hp_sum %>%
  pivot_wider(names_from = Time, values_from = c(mean, sd))

hp_sum <- hp_sum_wide %>%
  mutate(
    diff_After_Before = mean_After - mean_Before,
    percent_change = (diff_After_Before / mean_Before) * 100,
    diff_After_Before = sprintf("%+.0f", diff_After_Before),
    percent_change = sprintf("%+.0f%%", percent_change),
    Before = sprintf("%.0f (%.0f)", mean_Before, sd_Before),
    After = sprintf("%.0f (%.0f)", mean_After, sd_After)
  ) %>%
  select(Conditions, Before, After, diff_After_Before, percent_change)

hp_sum <- hp_sum %>%
  mutate(Conditions = if_else(str_detect(Conditions, "VLL"), 
                          paste0(Conditions, " (HP)"), 
                          Conditions))


# PL values
pl_var <- File_PL_data %>%
  select(-c(1:4)) %>%
  filter(minutes >= 17 & minutes <= 34) %>%
  mutate(Time = case_when(
    minutes >= 17 & minutes <= 19 ~ "Before",
    minutes >= 30 & minutes <= 34 ~ "After"))

pl_sum <- pl_var %>%
  group_by(Conditions, Time) %>%
  summarise(
    mean = mean(`F0`, na.rm = TRUE),
    sd = sd(`F0`, na.rm = TRUE),
    .groups = "drop")

pl_sum_wide <- pl_sum %>%
  pivot_wider(names_from = Time, values_from = c(mean, sd))

pl_sum <- pl_sum_wide %>%
  mutate(
    diff_After_Before = mean_After - mean_Before,
    percent_change = (diff_After_Before / mean_Before) * 100,
    diff_After_Before = sprintf("%+.0f", diff_After_Before),
    percent_change = sprintf("%+.0f%%", percent_change),
    Before = sprintf("%.0f (%.0f)", mean_Before, sd_Before),
    After = sprintf("%.0f (%.0f)", mean_After, sd_After)
  ) %>%
  select(Conditions, Before, After, diff_After_Before, percent_change)

pl_sum <- pl_sum %>%
  mutate(Conditions = if_else(str_detect(Conditions, "VLL"), 
                          paste0(Conditions, " (PL)"), 
                          Conditions))

# Post-stress resilience (60mins), slope

# HL 1 hour recovery
HL_slope <- File_HL_data[File_HL_data$minutes >= 43 & File_HL_data$minutes <= 105, ]

# Slope linear model
df_HL_slope <- HL_slope %>%
  group_by(Conditions) %>%
  do(modele = lm(`F0` ~ minutes, data = .)) %>%
  mutate(Slope = sprintf("%+.0f", coef(modele)[["minutes"]]))

# HP
HP_slope <- File_HP_data[File_HP_data$minutes >= 32 & File_HP_data$minutes <= 90, ]

df_HP_slope <- HP_slope %>%
  group_by(Conditions) %>%
  do(modele = lm(`F0` ~ minutes, data = .)) %>%
  mutate(Slope = sprintf("%+.0f", coef(modele)[["minutes"]]))

# PL
PL_slope <- File_PL_data[File_PL_data$minutes >= 30 & File_PL_data$minutes <= 93, ]

df_PL_slope <- PL_slope %>%
  group_by(Conditions) %>%
  do(modele = lm(`F0` ~ minutes, data = .)) %>%
  mutate(Slope = sprintf("%+.0f", coef(modele)[["minutes"]]))

# merge df_slope
df_slope <- bind_rows(df_HL_slope, df_HP_slope, df_PL_slope)

# Edit table
compiled_table <- bind_rows(hl_sum, hp_sum, pl_sum)
compiled_table$Slope <- df_slope$Slope # Add slope

compiled_table <- compiled_table %>%
  rename(
    `F0 After (a.u.)` = `After`,
    `F0 Before (a.u.)` = `Before`,
    `Δ After-Before (a.u.)` = `diff_After_Before`,
    `Changes (%)` = `percent_change`,
    `Recovery slope (F0.min-1)` = `Slope`)

compiled_table %>% 
  gt() %>%
  cols_width(
    everything() ~ px(200)
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = c(last_col(), last_col(offset = 1)))) 
