#Library packages ####
#install.packages("gt")
library(gt)
#install.packages("patchwork")
library(patchwork)
#install.packages("vegan")
library(vegan)
#install.packages("FactoMineR")
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
#install.packages("lattice")
library(lattice)
#install.packages("permute")
library(permute)
#install.packages("dplyr")
library(dplyr)
#install.packages("reshape2")
library(reshape2)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("ade4")
library(ade4)
#install.packages("readxl")
library(readxl)
#install.packages("plotly")
library(plotly)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("webshot2")
library(webshot2)
#install.packages("ggrepel")
library(ggrepel)

# FIGURE S1 ####

#ΣChla
df_spear_cor <- gdf_pigment %>%
  mutate('ΣChla' = rowSums(select(., starts_with("Chla"), starts_with("Pheopigments")))) %>%
  select(-starts_with("Chla"), -starts_with("Pheopigments"))
#Σpigments without Chla pigments
df_spear_cor <- df_spear_cor %>%
  mutate(`Σtot` = rowSums(select(., 5:31)))

correlation <- cor.test(df_spear_cor$ΣChla, df_spear_cor$`Σtot`, method = "spearman")
correlation_value <- correlation$estimate
R2_value <- correlation_value^2
p_value <- correlation$p.value


ggplot(df_spear_cor, aes(x = `Σtot`, y = ΣChla)) + 
  geom_point(color = "black") + 
  geom_smooth(method = "lm", se = FALSE, color = "#d62828") + 
  labs(
    x = "Total pigments (µg.gr sediment-1)", 
    y = "ΣChla (µg.gr sediment-1)", 
    title = "",
    subtitle = paste0(
      "Spearman correlation = ", round(correlation_value, 3), 
      ", p-value = ", format.pval(p_value, digits = 3), 
      ", R² = ", round(R2_value, 3)
    )
  ) + 
  theme_test()+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 14),
        axis.text.y = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15))
# FIGURE S2 ####

df_16 <- psmelt(final_biofilm16S)
df_16 <- df_16 %>%mutate(Phylum = ifelse(Phylum == "Proteobacteria", "Pseudomonadota", Phylum))

# Phylum group
df_agg_16 <- df_16 %>%
  group_by(Sample, Conditions, Biofilm_type, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# Abundance means
df_agg_16 <- df_agg_16 %>%
  group_by(Conditions, Phylum) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = 'drop')

# Relative abundances
df_agg_16 <- df_agg_16 %>%
  group_by(Phylum, Conditions) %>%
  summarise(Mean_Abundance = sum(Mean_Abundance), .groups = 'drop') %>%
  group_by(Conditions) %>%
  mutate(Relative_Abundance = Mean_Abundance / sum(Mean_Abundance)*100)

# Threshold
df_agg_16 <- df_agg_16 %>%
  mutate(Phylum = ifelse(Relative_Abundance < 1, "Others <1%", Phylum))

# Cleanning names
df_agg_16$Phylum <- gsub("\\s+", "", df_agg_16$Phylum)

# df subset
df_SB_16 <- subset(df_agg_16, grepl("^SB_", Conditions))
df_WSB_16 <- subset(df_agg_16, grepl("^WSB_", Conditions))

# Colors settings
all_Phylum_16 <- union(unique(df_SB_16$Phylum), unique(df_WSB_16$Phylum))
color_Phylum_16 <- colorRampPalette(c("black", "grey",   "#44aa99", "#01655E", "#c7eae5", "#c1a5cf","#dfc27d","steelblue", "#f3a482", "#4d4d4d", "#9970ab")) (14)
Phylum_colors_16 <- data.frame(Phylum = all_Phylum_16, Color = color_Phylum_16) 

# Phylum names
labels_italics <- setNames(
  lapply(all_Phylum_16, function(x) bquote(italic(.(x)))),
  all_Phylum_16)

# SB Plot
plot_SB_phylum_16 <- ggplot(df_SB_16, aes(x = Conditions, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_test() +
  labs(title = "SB", x = "", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = setNames(Phylum_colors_16$Color, Phylum_colors_16$Phylum),
                    labels = labels_italics)+
  theme(axis.title.y = element_text(size = 16, face = "bold"), 
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
plot_WSB_phylum_16 <- ggplot(df_WSB_16, aes(x = Conditions, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_test() +
  labs(title = "WSB", x = "", y = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = setNames(Phylum_colors_16$Color, Phylum_colors_16$Phylum))+
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

(plot_SB_phylum_16 + plot_WSB_phylum_16)+
  plot_annotation(
    title = "",
    caption = "Stress Conditions") & theme(plot.caption = element_text(hjust = 0.5, size = 16, face = "bold"),
                                           plot.title = element_text(face="bold", size="20"))

# Average pigment composition across all samples
barplot_16_all <- rbind(df_SB_16, df_WSB_16)
aggregate(Relative_Abundance ~ Phylum, data = barplot_16_all, mean)


# FIGURE S3 ####

# df structure
metadata_mfa <- df_community[, 1:3]             
taxonomy18s_mfa <- df_community[, 4:3355]  
taxonomy16s_mfa <- df_community[, 3356:10657] 

# df
mfa_data <- df_community[, 4:10657]

# factor
Groups <- as.factor(paste(metadata_mfa$Groups)) 

mfa_result <- MFA(mfa_data, 
                  group = c(ncol(taxonomy18s_mfa), ncol(taxonomy16s_mfa)), 
                  type = c("s","s"), # standardized qualitative data 's'
                  name.group = c("Eukaryote", "Procaryote"))

# ConvexHull function
source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R")

# Plot
plot_mfa <- fviz_mfa(mfa_result, "all",
         axes = c(1, 2),
         habillage = Groups,
         label = "none") +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(), linetype = "none",  fill = "none") +
  labs(color = "Groups") +
  scale_fill_manual(values = c("SB" = "#CC6677", "WSB" = "#81B29A", "Eukaryote" = "#073B4C", "Procaryote"="orange")) +
  scale_color_manual(values = c("SB" = "#CC6677", "WSB" = "#81B29A", "Eukaryote" = "#073B4C", "Procaryote"="orange")) +
  geom_convexhull(aes(fill = Groups), alpha = 0.3) +
  theme_test()+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 14),
        axis.text.y = element_text(color = "black", face = "bold", size = 14),
        legend.title = element_text(face = "bold"))

ggplotly(plot_mfa)

# Anosim between SB and WSB community
coord_axes <- mfa_result$ind$coord[,1]
groups_fac <- df_community$Groups
coord_axes <- data.frame(coord_axes, group = groups_fac)

# Euclidean matrix
distance_matrix <- dist(coord_axes[,1], method = "euclidean")
anosim(distance_matrix, coord_axes$group, permutations = 9999)

# FIGURE S4 ####

  # PROKARYOTES

richesse_spe16 <- estimate_richness(biofilm16S)
richesse_spe16$Pielou <- richesse_spe16$Shannon / log(richesse_spe16$Observed) # Pielou

# Metadata
metadata_div<- df_community[, 1:3]  
rownames(metadata_div) <- metadata_div$File
richesse_spe16 <- merge(richesse_spe16, metadata_div, by = "row.names", all.x = TRUE)
richesse_spe16 <- richesse_spe16 %>%mutate(Treatment = paste(Conditions, Groups, sep = "_"))

# Df format
richesse_spe16_long <- richesse_spe16 %>%
  select(Observed, Chao1, Shannon, InvSimpson, Pielou, Groups, Treatment) %>%
  gather(key = "Indice", value = "Valeur", Observed, Chao1, Shannon, InvSimpson, Pielou)

# Order
richesse_spe16_long$Indice <- factor(richesse_spe16_long$Indice, levels = c("Observed", "Chao1", "Shannon", "InvSimpson", "Pielou"))

# Colors settings
colors_groups <- c("SB" = "#CC6677", "WSB" = "#81B29A")

# Stat comparisons 
my_comparisons <- list(c("WSB", "SB"))

# Plot
plot_16_div <- ggplot(richesse_spe16_long, aes(x = Groups, y = Valeur, fill = Groups, color = Groups)) +
  geom_boxplot(outlier.shape = NA, size = 0.5, alpha = 1) + 
  geom_jitter(width = 0.2, shape = 16, size = 2, alpha = 1) + 
  scale_fill_manual(values = colors_groups) + 
  scale_color_manual(values = c("SB" = "black", "WSB" = "black")) + 
  facet_wrap(~ Indice, scales = "free_y") +  
  labs(x = "", y = "Index values", title = "Prokaryotic community") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif")+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 14),
        axis.text.y = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))

  # EUKARYOTES

richesse_spe18 <- estimate_richness(biofilm18S)
richesse_spe18$Pielou <- richesse_spe18$Shannon / log(richesse_spe18$Observed) # Pielou

# Metadata
richesse_spe18 <- merge(richesse_spe18, metadata_div, by = "row.names", all.x = TRUE)
richesse_spe18 <- richesse_spe18 %>%mutate(Treatment = paste(Conditions, Groups, sep = "_"))

# Df format
richesse_spe18_long <- richesse_spe18 %>%
  select(Observed, Chao1, Shannon, InvSimpson, Pielou, Groups, Treatment) %>%
  gather(key = "Indice", value = "Valeur", Observed, Chao1, Shannon, InvSimpson, Pielou)

# Order
richesse_spe18_long$Indice <- factor(richesse_spe18_long$Indice, levels = c("Observed", "Chao1", "Shannon", "InvSimpson", "Pielou"))

# Plot
plot_18_div <- ggplot(richesse_spe18_long, aes(x = Groups, y = Valeur, fill = Groups, color = Groups)) +
  geom_boxplot(outlier.shape = NA, size = 0.5, alpha = 1) + 
  geom_jitter(width = 0.2, shape = 16, size = 2, alpha = 1) + 
  scale_fill_manual(values = colors_groups) + 
  scale_color_manual(values = c("SB" = "black", "WSB" = "black")) + 
  facet_wrap(~ Indice, scales = "free_y") +  
  labs(x = "", y = "Index values", title = "Eukaryotic community") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif")+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 14),
        axis.text.y = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))

(plot_18_div + plot_16_div) + 
  plot_annotation(
    title = "",
    caption = "Groups") & theme(plot.caption = element_text(hjust = 0.5, size = 20, face = "bold"),
                                    plot.title = element_text(size="20"))

# FIGURE S5 ####

  # (A)

# 1/2 day informations
highlight_times <- Engogenous_mig_df[Engogenous_mig_df$hour %in% c("00:00:00", "12:00:00"), ]

# Plot
a <- ggplot() +
  # Sealevel
  geom_line(data = Engogenous_mig_df, aes(x = Time, y = SeaLevel * 1000), color = "#5fa8d3", size = 0.8) +  
  scale_y_continuous(
    name = "Surface microalgal biomass (F0 a.u.)", 
    sec.axis = sec_axis(
      trans = ~ . / 1000, # scale
      name = "In Situ sea level (m)",  
      breaks = seq(0, 10, by = 1),
      labels = seq(0, 10, by = 1)   
    )
  ) +
  # FO
  geom_line(data = Engogenous_mig_df, aes(x = Time, y = `F0`), color = "black") +
  geom_point(data = Engogenous_mig_df, aes(x = Time, y = `F0`, color = Daylight), size = 2) +
  scale_color_manual(values = c("black", "darkgrey"), labels = c("Night", "Day")) +
  labs(x = "Time (Days)", y = "", 
       title = "A",
       color = "Photoperiod") +
  theme_test() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size = 18), 
    axis.title.y.right = element_text(color = "#5fa8d3", size = 16, face = "bold"),
    axis.text.y.right = element_text(color = "#5fa8d3", size = 14, face = "bold"),   
    panel.grid.major = element_line(colour = "white"),  
    panel.grid.minor = element_line(colour = "white"),
    legend.position = "bottom") +
  geom_vline(data = highlight_times, aes(xintercept = Time), linetype = "dashed", color = "grey")

  # (B)

#ACF
rownames(Engogenous_mig_df) <- seq_len(nrow(Engogenous_mig_df))
x<-Engogenous_mig_df$`F0`
result_acf<-acf(x,lag.max=1500,type = "correlation", plot = F)
acf_data <- data.frame(Lag = result_acf$lag, ACF = result_acf$acf)

# Plot
b <- ggplot(acf_data, aes(x = Lag, y = ACF)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_segment(aes(xend = Lag, yend = 0), color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(0.05, -0.05), linetype = "dashed", color = "#d62828") +
  labs(
    x = "Lag",
    y = "Autocorrelation",
    title = "B") +
  theme_test() +
  theme(
    axis.text.x = element_text(color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", size = 18))+
  annotate("text", 
           x = 534, y = 1,
           label = "Average maximum periodicity: 47", 
           color = "black", 
           size = 5, 
           fontface = "bold")+
  annotate("text", 
           x = 540, y = 0.85,
           label = "Measurement every 30 minutes", 
           color = "black", 
           size = 5, 
           fontface = "bold")+
  annotate("text", 
           x = 525, y = 0.70, 
           label = "Average migration periodicity: 23.6h", 
           color = "black", 
           size = 5, 
           fontface = "bold")


# Period
print(data.frame(result_acf$lag,result_acf$acf))

p1<-47*30/60 # lag dist between 2 max, 30 because lag of `F0` mesures, 60 minutes for 1h 
p2<-48*30/60 
p3<-47*30/60
p4<-47*30/60
p5<-47*30/60

# Mean period
(p1+p2+p3+p4+p5)/5

  # (C)

  # wilcox.test `F0` ~ Daylight

# Cut 1/2 day start and end time serie
cor_DL_data <- Engogenous_mig_df %>%slice((17):(n() - 17))

cor_DL_data$Daylight<-gsub("TRUE","1",cor_DL_data$Daylight)
cor_DL_data$Daylight<-gsub("FALSE","0",cor_DL_data$Daylight)
cor_DL_data$Daylight<-as.numeric(cor_DL_data$Daylight)

wilcox.test(`F0` ~ Daylight, alternative = "less", data = cor_DL_data) 

  # Spearman correlation, average `F0` ~ average SeaLevel

# Cut 1/2 day start and end time serie
SL_data <- Engogenous_mig_df %>%slice((17):(n() - 17))

# Days factor
d<-factor(substr(SL_data$Time,1,10),levels=c("2023-12-23","2023-12-24","2023-12-25","2023-12-26","2023-12-27","2023-12-28","2023-12-29","2023-12-30","2023-12-31","2024-01-01","2024-01-02","2024-01-03"))

# F0 mean by days
av_F0<-aggregate(SL_data$`F0`,by=list(d),mean)
av.frame.F0<-data.frame(TIME.av.F0=as.POSIXct(strptime(paste(av_F0$Group.1,"12:00:00"),format="%Y-%m-%d %H:%M:%S")),F0.av=av_F0$x)

# SL mean by days
av_SL<-aggregate(SL_data$SeaLevel,by=list(d),mean)
av.frame.SL<-data.frame(TIME.av.SL=as.POSIXct(strptime(paste(av_SL$Group.1,"12:00:00"),format="%Y-%m-%d %H:%M:%S")),SL.av=av_SL$x)

# Spearman cor
cor_SL_data <- av.frame.F0 %>%
  rename(Time = TIME.av.F0) %>%
  inner_join(av.frame.SL %>% rename(Time = TIME.av.SL), by = "Time")

cor_SL <- cor.test(av.frame.SL$SL.av, av.frame.F0$F0.av, method = 'spearman')

# Retrieve data
cor_SL_values <- round(cor_SL$estimate, 2)
R2_value_SL <- round(cor_SL_values^2, 2)
p_value_SL <- cor_SL$p.value

# Plot
c <- ggplot(cor_SL_data, aes(x = SL.av, y = F0.av)) +
  geom_point(color = "black") + 
  geom_smooth(method = "lm", se = FALSE, color = "#d62828") +
  labs(
    x = "Daily average sea level (meters)", 
    y = "Daily average biomass (F0 a.u.)", 
    title = "C",
    subtitle = paste0(
      "Spearman correlation = ", round(cor_SL_values, 3), 
      ", p-value = ", format.pval(p_value_SL, digits = 3), 
      ", R² = ", round(R2_value_SL, 3))) + 
  theme_test()+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 14),
        axis.text.y = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.title = element_text(face = "bold", size = 18))

a / b / c + 
plot_layout(heights = c(2, 1, 1))

# FIGURE S6 ####

df_barplot_pigment <- gdf_pigment_perc

#Carotenoid-like group
df_barplot_pigment <- df_barplot_pigment %>%
  mutate(`Unknown-carotenoids` = rowSums(select(., starts_with("Uc")))) %>%
  select(-starts_with("Uc"))

#fucoxanthin-like group
df_barplot_pigment<-df_barplot_pigment %>%
  mutate(`Fucoxanthin-likes` = rowSums(select(., starts_with("Fl")))) %>%
  select(-starts_with("Fl"))

colors_barplot_pigments <- c(
  "V" = "#774993", 
  "A" = "#1f77b4", 
  "Z" = "#93375A", 
  "Dt" = "#0A0701", 
  "Dd" = "#c5b0d5",
  "Chla" = "#4C8543", 
  "Chla-derivatives" = "#B3EE98", 
  "Pheopigments" = "darkgreen",
  "F" = "#81724F", 
  "Fucoxanthin-likes" = "#c49c94", 
  "Unknown-carotenoids" = "#ff9896", 
  "Chlc2" = "#ffbb78", 
  "ββ" = "#ff7f0e", 
  "βε" = "#aec7e8", 
  "L" = "#dbdb8d", 
  "N" = "#382469")

# Subset by groups
df_barplot_pigment_sb <- subset(df_barplot_pigment, Groups == "SB")
df_barplot_pigment_wsb <- subset(df_barplot_pigment, Groups == "WSB")

#SB
data_perc_sb <- df_barplot_pigment_sb %>%
  select(Conditions, 5:20) %>%
  pivot_longer(cols = -Conditions, names_to = "Pigments", values_to = "Percentage")

# Plot
data_perc_sb$Conditions <- factor(data_perc_sb$Conditions, levels = c("VLL", "HL", "HP", "PL"))

barplot_pigment1 <- ggplot(data_perc_sb, aes(x = Conditions, y = Percentage, fill = Pigments)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_test() +
  labs(title = "SB", x = "", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = colors_barplot_pigments)+
  theme(axis.title.y = element_text(size = 20, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "bottom")+
  annotate("text", 
           x = 1, y = -0.05, 
           label = "VLL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 2, y = -0.05, 
           label = "HL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 3, y = -0.05, 
           label = "HP", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 4, y = -0.05, 
           label = "PL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")

#WSB
data_perc_wsb <- df_barplot_pigment_wsb %>%
  select(Conditions, 5:20) %>%
  pivot_longer(cols = -Conditions, names_to = "Pigments", values_to = "Percentage")

# Plot
data_perc_wsb$Conditions <- factor(data_perc_wsb$Conditions, levels = c("VLL", "HL", "HP", "PL"))

barplot_pigment2 <- ggplot(data_perc_wsb, aes(x = Conditions, y = Percentage, fill = Pigments)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_test() +
  labs(title = "WSB", x = "", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = colors_barplot_pigments)+
  theme(axis.title.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none")+
  annotate("text", 
           x = 1, y = -0.05, 
           label = "VLL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 2, y = -0.05, 
           label = "HL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 3, y = -0.05, 
           label = "HP", 
           size = 6, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 4, y = -0.05, 
           label = "PL", 
           size = 6, 
           fontface = "bold", 
           colour = "black")


(barplot_pigment1 + barplot_pigment2) + 
  plot_annotation(
    title = "",
    caption = "Conditions") & theme(plot.caption = element_text(hjust = 0.5, size = 20, face = "bold"),
                                    plot.title = element_text(size="20"))

# Average pigment composition across all samples
barplot_pig_all <- rbind(data_perc_sb, data_perc_wsb)
aggregate(Percentage ~ Pigments, data = barplot_pig_all, mean)

# FIGURE S7 ####

bca_sb <- subset(df_pigment_perc, Groups == "SB")
bca_wsb <- subset(df_pigment_perc, Groups == "WSB")

source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R") # Convex function

#SB BCA
res.pca_sb <- dudi.pca(bca_sb [, - c(1:4)], scannf = FALSE, nf = 2) #PCA

res.bca.pigment_sb <- bca(res.pca_sb, as.factor(bca_sb$Conditions), scannf = FALSE, nf = 2)
res.bca.pigment_sb$ratio 

# eig values
varexp2<-res.bca.pigment_sb$eig*100/sum(res.bca.pigment_sb$eig)

fac<-as.factor(bca_sb$Conditions)

# Plot
bca1 <- ggplot(res.bca.pigment_sb$ls,aes(x=res.bca.pigment_sb$ls[,1],y=res.bca.pigment_sb$ls[,2],col=fac))+
  geom_point()+
  scale_fill_manual(values=c("VLL" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2")) +
  scale_color_manual(values=c("VLL" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2")) +
  geom_convexhull(alpha = 0.3,aes(fill = fac))+
  xlab(paste("Axis ",1," : ",round(varexp2[1],2),"%"))+
  ylab(paste("Axis ",2," : ",round(varexp2[2],2),"%"))+
  labs(fill="Explanatory factors",col="Explanatory factors")+
  theme_test() +
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 16),
        axis.text.y = element_text(color = "black", face = "bold", size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.position = "bottom")+
  annotate("text", 
           x = -5.6, y = 5, 
           label = "SB (84%)", 
           size = 8, 
           fontface = "bold", 
           colour = "black")


#WSB BCA
res.pca.ade4 <- dudi.pca(bca_wsb [, - c(1:4)], scannf = FALSE, nf = 2) # PCA
res.bca.pigment_wsb <- bca(res.pca.ade4, as.factor(bca_wsb$Conditions), scannf = FALSE, nf = 2)
res.bca.pigment_wsb$ratio 

varexp2<-res.bca.pigment_wsb$eig*100/sum(res.bca.pigment_wsb$eig)

fac<-as.factor(bca_wsb$Conditions)

# Plot
bca2 <- ggplot(res.bca.pigment_wsb$ls,aes(x=res.bca.pigment_wsb$ls[,1],y=res.bca.pigment_wsb$ls[,2],col=fac))+
  geom_point()+
  scale_fill_manual(values=c("VLL" = "black", "HP" = "steelblue", "HL" = "#E7B800", "PL" = "#B869B8")) +
  scale_color_manual(values=c("VLL" = "black", "HP" = "steelblue", "HL" = "#E7B800", "PL" = "#B869B8")) +
  geom_convexhull(alpha = 0.3,aes(fill = fac))+
  xlab(paste("Axis ",1," : ",round(varexp2[1],2),"%"))+
  ylab(paste("Axis ",2," : ",round(varexp2[2],2),"%"))+
  labs(fill="Explanatory factors",col="Explanatory factors")+
  theme_test() +
  
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 16),
        axis.text.y = element_text(color = "black", face = "bold", size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.position = "none")+
  annotate("text", 
           x = -3.5, y = 5, 
           label = "WSB (66%)", 
           size = 8, 
           fontface = "bold", 
           colour = "black")

(bca1 + bca2) + 
  plot_annotation(
    title = "",
    caption = "Conditions") & theme(plot.caption = element_text(hjust = 0.5, size = 20, face = "bold"),
                                    plot.title = element_text(size="20"))

  # Anosim, impact of oxidative conditions

# SB samples
df_pigment_perc_sb <- subset(df_pigment_perc, Groups == "SB")

pigments_anosim_sb <- df_pigment_perc_sb %>% select(5:ncol(df_pigment_perc_sb))
pca_pigments_anosim_sb <- prcomp(pigments_anosim_sb, scale = TRUE) # PCA

coord_axes <- pca_pigments_anosim_sb$x[, 1:2]
group_fac <- df_pigment_perc_sb$Conditions

coord_axes <- data.frame(coord_axes, group = group_fac)

# Euclidean distance matrix
distance_matrix_pigment <- dist(coord_axes[, c("PC1", "PC2")], method = "euclidean")
anosim(distance_matrix_pigment, coord_axes$group, permutations = 9999)


# WSB samples
df_pigment_perc_wsb <- subset(df_pigment_perc, Groups == "WSB")

pigments_anosim_wsb <- df_pigment_perc_wsb %>% select(5:ncol(df_pigment_perc_wsb))
pca_pigments_anosim_wsb <- prcomp(pigments_anosim_wsb, scale = TRUE) # PCA

coord_axes <- pca_pigments_anosim_wsb$x[, 1:2]
group_fac <- df_pigment_perc_wsb$Conditions

coord_axes <- data.frame(coord_axes, group = group_fac)

# Euclidean distance matrix
distance_matrix_pigment <- dist(coord_axes[, c("PC1", "PC2")], method = "euclidean")
anosim(distance_matrix_pigment, coord_axes$group, permutations = 9999)

# FIGURE S8 ####

# SB

  #Table
df_pigment_perc_sb <- subset(df_pigment_perc, Groups == "SB")

simper_sb <- df_pigment_perc_sb[, -c(1:4)]
simper_sb_var <- df_pigment_perc_sb[, -c(5:48)]
sim_sb <- simper(simper_sb, simper_sb_var$Conditions, permutations = 9999)

sim_sb_df <- bind_rows(
  lapply(names(sim_sb), function(name) {
    sim <- sim_sb[[name]]
    df <- as.data.frame(sim)
    df$species <- rownames(df)
    df <- df[, c("species", "average", "sd", "p")]
    df$contrast <- name
    df[, c("contrast", "species", "average", "sd", "p")]
  })
)

sim_sb_df <- sim_sb_df %>%filter(p > 1e-7 & p < 0.05) # Only significant values

sim_sb_df %>%
  gt() %>%
  tab_header(
    title = "SIMPER SB results - Stress Contrast"
  ) %>%
  cols_label(
    contrast = "Contrast",
    species = "Species",
    average = "Mean difference",
    sd = "sd",
    p = "p-value"
  ) %>%
  fmt_number(
    columns = c(average, sd, p),
    decimals = 4
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()))%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))

# Pigment number impacted
sim_sb_df_contrast <- sim_sb_df %>%
  count(contrast)

  # Plot

pigments_sb <- df_pigment_perc_sb %>% select(5:ncol(df_pigment_perc_sb))
pca_pigments_sb <- prcomp(pigments_sb, scale = TRUE) # PCA

data_sb <- cbind(df_pigment_perc_sb, pca_pigments_sb$x) # Add PCA results

n<-45 # Scale

source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R") # Convex function


biplot_sb <- ggplot(data_sb, aes(x = -PC1, y = -PC2, color = Conditions)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("VLL" = "#324851", "HP" = "#5C88C4", "HL"="#E7B800", "PL"="#937DC2")) +
  scale_fill_manual(values = c("VLL" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2"))+
  geom_convexhull(aes(fill = Conditions), alpha=0.6)+
  annotate("segment",
           x = 0, y = 0,
           xend = -pca_pigments_sb$rotation[,1] * n,
           yend = -pca_pigments_sb$rotation[,2] * n,
           arrow = arrow(length = unit(0.2, "cm")),
           color = "darkgrey",
           linewidth = 0.2) +
  geom_text_repel(
    data = as.data.frame(pca_pigments_sb$rotation),
    aes(x = -PC1 * n, y = -PC2 * n, label = colnames(pigments_sb)),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = 100) +
  labs(
    title = "A",
    x = paste0("PC 1 (", round(100 * pca_pigments_sb$sdev[1]^2 / sum(pca_pigments_sb$sdev^2), 1), "%)"),
    y = paste0("PC 2 (", round(100 * pca_pigments_sb$sdev[2]^2 / sum(pca_pigments_sb$sdev^2), 1), "%)")) +
  guides(
    color = guide_legend(override.aes = list(shape = NA)),
    shape = guide_legend(),
    fill = guide_legend(override.aes = list(shape = c(19, 17), color = "black", linetype = "dashed"))) +
  theme_test() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none")


# WSB

#Table
df_pigment_perc_wsb <- subset(df_pigment_perc, Groups == "WSB")

simper_wsb <- df_pigment_perc_wsb[, -c(1:4)]
simper_wsb_var <- df_pigment_perc_wsb[, -c(5:48)]
sim_wsb <- simper(simper_wsb, simper_wsb_var$Conditions, permutations = 9999)

sim_wsb_df <- bind_rows(
  lapply(names(sim_wsb), function(name) {
    sim <- sim_wsb[[name]]
    df <- as.data.frame(sim)
    df$species <- rownames(df)
    df <- df[, c("species", "average", "sd", "p")]
    df$contrast <- name
    df[, c("contrast", "species", "average", "sd", "p")]
  })
)

sim_wsb_df <- sim_wsb_df %>%filter(p > 1e-7 & p < 0.05) # Only significant values

sim_wsb_df %>%
  gt() %>%
  tab_header(
    title = "SIMPER WSB results - Stress Contrast"
  ) %>%
  cols_label(
    contrast = "Contrast",
    species = "Species",
    average = "Mean difference",
    sd = "sd",
    p = "p-value"
  ) %>%
  fmt_number(
    columns = c(average, sd, p),
    decimals = 4
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()))%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))

# Pigment number impacted
sim_wsb_df_contrast <-sim_wsb_df %>%
  count(contrast)

# Plot

pigments_wsb <- df_pigment_perc_wsb %>% select(5:ncol(df_pigment_perc_wsb))
pca_pigments_wsb <- prcomp(pigments_wsb, scale = TRUE) # PCA

data_wsb <- cbind(df_pigment_perc_wsb, pca_pigments_wsb$x) # Add PCA results

n<-25 # Scale

biplot_wsb <- ggplot(data_wsb, aes(x = -PC1, y = -PC2, color = Conditions)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("VLL" = "#324851", "HP" = "#5C88C4", "HL"="#E7B800", "PL"="#937DC2")) +
  scale_fill_manual(values = c("VLL" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2"))+
  geom_convexhull(aes(fill = Conditions), alpha=0.6)+
  annotate("segment",
           x = 0, y = 0,
           xend = -pca_pigments_wsb$rotation[,1] * n,
           yend = -pca_pigments_wsb$rotation[,2] * n,
           arrow = arrow(length = unit(0.2, "cm")),
           color = "darkgrey",
           linewidth = 0.2) +
  geom_text_repel(
    data = as.data.frame(pca_pigments_wsb$rotation),
    aes(x = -PC1 * n, y = -PC2 * n, label = colnames(pigments_wsb)),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = 100) +
  labs(
    title = "B",
    x = paste0("PC 1 (", round(100 * pca_pigments_wsb$sdev[1]^2 / sum(pca_pigments_wsb$sdev^2), 1), "%)"),
    y = paste0("PC 2 (", round(100 * pca_pigments_wsb$sdev[2]^2 / sum(pca_pigments_wsb$sdev^2), 1), "%)")) +
  guides(
    color = guide_legend(override.aes = list(shape = NA)),
    shape = guide_legend(),
    fill = guide_legend(override.aes = list(shape = c(19, 17), color = "black", linetype = "dashed"))) +
  theme_test() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right")

# Final graph
sim_wsb_df_contrast$Group <- "WSB"
sim_sb_df_contrast$Group <- "SB"
values_df <- rbind(sim_sb_df_contrast, sim_wsb_df_contrast)
values_df <- values_df[values_df$contrast %in% c("HL_VLL", "HP_VLL", "PL_VLL"), ]
values_df$Percentage <- (values_df$n / 44) * 100
values_df <- values_df %>%rename(Contrast = contrast)

pig_impact <- ggplot(values_df, aes(x = Contrast, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("SB" = "#CC6677", "WSB" = "#81B29A")) +
  theme_test(base_size = 16) +
  labs(
    title = "C",
    x = "Contrast",
    y = "Percentage",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right")


biplot_sb / biplot_wsb / pig_impact + 
  plot_layout(heights = c(2, 2, 1))

# FIGURE S9 ####

# Distance matrix
coord_axes <- pca_pigments$x[,1:2]  # Run FIGURE 5 before !
group_fac <- df_pigment_perc$Groups # metadata
coord_axes <- data.frame(coord_axes, group = group_fac)
distance_matrix_pigment <- dist(coord_axes[, 1:2], method = "euclidean") # Euclidean distance matrix

distance_matrix_pigment = as.matrix(distance_matrix_pigment)
treatment_fac = split(seq_along(df_pigment_perc$Treatment), df_pigment_perc$Treatment)

comparisons <- list(
  "HL" = distance_matrix_pigment[treatment_fac$VLL_SB, treatment_fac$HL_SB],
  "PL" = distance_matrix_pigment[treatment_fac$VLL_SB, treatment_fac$PL_SB],
  "HP" = distance_matrix_pigment[treatment_fac$VLL_SB, treatment_fac$HP_SB],
  "HL " = distance_matrix_pigment[treatment_fac$VLL_WSB, treatment_fac$HL_WSB],
  "PL " = distance_matrix_pigment[treatment_fac$VLL_WSB, treatment_fac$PL_WSB],
  "HP " = distance_matrix_pigment[treatment_fac$VLL_WSB, treatment_fac$HP_WSB])

comparison_stats <- lapply(comparisons, function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x), sd = sd(x))
})

# Display order
stress_order <- c("HL", "HP", "PL", "HL ", "HP ", "PL ")
comparison_stats_ordered <- comparison_stats[stress_order]

# Mean & sd values
means <- sapply(comparison_stats_ordered, function(x) x["mean"])
sds <- sapply(comparison_stats_ordered, function(x) x["sd"])

# df
table_ave_dist_pigment<- data.frame(Comparaison = names(means), Mean = means, SD = sds)
rownames(table_ave_dist_pigment) <- NULL
table_ave_dist_pigment$Comparaison <- sub("\\.mean$", "", table_ave_dist_pigment$Comparaison)

# xlab order
table_ave_dist_pigment$Comparaison <- factor(table_ave_dist_pigment$Comparaison, levels = stress_order)

# Plot
ggplot(table_ave_dist_pigment, aes(x = Comparaison, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                width = 0.2, color = "darkgrey", size = 0.7) +
  geom_point() +
  geom_text(aes(label = round(Mean, 2)), 
            vjust = -1, color = "black", size = 5) + 
  labs(
    title = "",
    x = "Conditions",
    y = "Average Euclidean dissimilarities (± sd)"
  ) +
  theme_test() +
  
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 16),
        axis.text.y = element_text(color = "black", face = "bold", size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"))+
  
  annotate("segment", 
           x = 0.51, xend = 3.48,  
           y = 7.5, yend = 7.5,  
           colour = "black", 
           size = 1) +
  
  annotate("segment", 
           x = 3.51, xend = 6.5,  
           y = 7.5, yend = 7.5,  
           colour = "black", 
           size = 1) +
  annotate("text", 
           x = 2, y = 7.8, 
           label = "SB", 
           size = 8, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 5, y = 7.8, 
           label = "WSB", 
           size = 8, 
           fontface = "bold", 
           colour = "black")


# FIGURE S10 ####

df_heatmap_pigment <- gdf_pigment

df_heatmap_pigment <- df_heatmap_pigment %>%mutate('ΣChla' = rowSums(select(., starts_with("Chla"), starts_with("Pheo")))) %>%select(-starts_with("Chla"), -starts_with("Pheo"))

df_heatmap_pigment[, 5:32] <- sweep(df_heatmap_pigment[, 5:32],MARGIN = 1,STATS = df_heatmap_pigment$`ΣChla`,FUN = "/") # pigment concentration (µg pigments.g of dry matter-1) by chla unit conversion

# Subset
df_heatmap_pigment_sb <- subset(df_heatmap_pigment, grepl("_SB_", File))
df_heatmap_pigment_wsb <- subset(df_heatmap_pigment, grepl("_WSB_", File))

SB_VLL <- df_heatmap_pigment_sb %>% filter(Conditions == "VLL")
SB_HL <- df_heatmap_pigment_sb %>% filter(Conditions == "HL")
SB_PL <- df_heatmap_pigment_sb %>% filter(Conditions == "PL")
SB_HP <- df_heatmap_pigment_sb %>% filter(Conditions == "HP")

WSB_VLL <- df_heatmap_pigment_wsb %>% filter(Conditions == "VLL")
WSB_HL <- df_heatmap_pigment_wsb %>% filter(Conditions == "HL")
WSB_PL <- df_heatmap_pigment_wsb %>% filter(Conditions == "PL")
WSB_HP <- df_heatmap_pigment_wsb %>% filter(Conditions == "HP")

# wilcox.test loop function
perform_wilcox_pv_SB <- function(group1, group2) {
  p_values <- c()
  for (param in colnames(df_heatmap_pigment_sb)[5:31]) {
    wilcox_test <- wilcox.test(group1[[param]], group2[[param]])
    p_values <- c(p_values, wilcox_test$p.value)
  }
  return(p_values)
}

perform_wilcox_test_pv_WSB <- function(group1, group2) {
  p_values <- c()
  for (param in colnames(df_heatmap_pigment_wsb)[5:31]) {
    wilcox_test <- wilcox.test(group1[[param]], group2[[param]])
    p_values <- c(p_values, wilcox_test$p.value)
  }
  return(p_values)
}

# p-values
p_values_pigment_HL_SB <- perform_wilcox_pv_SB(SB_VLL, SB_HL)
p_values_pigment_PL_SB <- perform_wilcox_pv_SB(SB_VLL, SB_PL)
p_values_pigment_HP_SB <- perform_wilcox_pv_SB(SB_VLL, SB_HP)

p_values_pigment_HL_WSB <- perform_wilcox_test_pv_WSB(WSB_VLL, WSB_HL)
p_values_pigment_PL_WSB <- perform_wilcox_test_pv_WSB(WSB_VLL, WSB_PL)
p_values_pigment_HP_WSB <- perform_wilcox_test_pv_WSB(WSB_VLL, WSB_HP)

# p-values df
df_heatmap_pigments_pv <- data.frame(
  "HL" = p_values_pigment_HL_SB,
  "PL" = p_values_pigment_PL_SB,
  "HP" = p_values_pigment_HP_SB,
  "HL " = p_values_pigment_HL_WSB,
  "PL " = p_values_pigment_PL_WSB,
  "HP " = p_values_pigment_HP_WSB)

# colnames
df_heatmap_pigments_pv <- t(df_heatmap_pigments_pv)
colnames(df_heatmap_pigments_pv) <- colnames(df_heatmap_pigment_wsb)[5:31]

# df format
df_heatmap_pigment_pv_long <- melt(as.matrix(df_heatmap_pigments_pv))
colnames(df_heatmap_pigment_pv_long) <- c("Condition", "Pigment", "Value")

# p-values categories
df_heatmap_pigment_pv_long$Categorie <- cut(df_heatmap_pigment_pv_long$Value,breaks = c(-Inf, 0.05, 0.1, Inf),labels = c("< 0.05", "0.05 < > 0.1", "> 0.1"),right = FALSE)

# Mean values loop function
perform_mean_pigment_sb <- function(group1, group2) {
  mean_diff <- c()
  for (param in colnames(df_heatmap_pigment_sb)[5:31]) {
    diff <- mean(group2[[param]]) - mean(group1[[param]])
    mean_diff <- c(mean_diff, diff)
  }
  return(mean_diff)
}

perform_mean_pigment_wsb <- function(group1, group2) {
  mean_diff <- c()
  for (param in colnames(df_heatmap_pigment_wsb)[5:31]) {
    diff <- mean(group2[[param]]) - mean(group1[[param]])
    mean_diff <- c(mean_diff, diff)
  }
  return(mean_diff)
}

# Means calculations stress vs VLL
mean_diff_pigment_HL_SB <- perform_mean_pigment_sb(SB_VLL, SB_HL)
mean_diff_pigment_PL_SB <- perform_mean_pigment_sb(SB_VLL, SB_PL)
mean_diff_pigment_HP_SB <- perform_mean_pigment_sb(SB_VLL, SB_HP)

mean_diff_pigment_HL_WSB <- perform_mean_pigment_wsb(WSB_VLL, WSB_HL)
mean_diff_pigment_PL_WSB <- perform_mean_pigment_wsb(WSB_VLL, WSB_PL)
mean_diff_pigment_HP_WSB <- perform_mean_pigment_wsb(WSB_VLL, WSB_HP)

# Means df
df_heatmap_pigment_mean_diff <- data.frame(
  "HL" = mean_diff_pigment_HL_SB,
  "PL" = mean_diff_pigment_PL_SB,
  "HP" = mean_diff_pigment_HP_SB,
  "HL " = mean_diff_pigment_HL_WSB,
  "PL " = mean_diff_pigment_PL_WSB,
  "HP " = mean_diff_pigment_HP_WSB)

# Colnames
df_heatmap_pigment_mean_diff <- t(df_heatmap_pigment_mean_diff)
colnames(df_heatmap_pigment_mean_diff) <- colnames(df_heatmap_pigment_sb)[5:31]

# Df format
df_heatmap_pigment_mean_diff_long <- melt(as.matrix(df_heatmap_pigment_mean_diff))
colnames(df_heatmap_pigment_mean_diff_long) <- c("Condition", "Pigment", "Différence")

# Add arrows
df_heatmap_pigment_mean_diff_long$Arrow <- ifelse(df_heatmap_pigment_mean_diff_long$`Différence` > 0, "↑", "↓")

# Add arrow colors
df_heatmap_pigment_mean_diff_long$Arrow_Color <- ifelse(df_heatmap_pigment_mean_diff_long$Arrow == "↓", "#d62828", "black")
df_heatmap_pigment_pv_long <- cbind(df_heatmap_pigment_pv_long, df_heatmap_pigment_mean_diff_long[c("Arrow", "Arrow_Color")]) # merge

# Plot
ggplot(df_heatmap_pigment_pv_long, aes(x = Condition, y = Pigment, fill = Categorie)) + 
  geom_tile(color = "black", size = 0.5) + 
  scale_fill_manual(name = "P-values Mann-Whitney U test", 
                    values = c("< 0.05" ="#6096ba", "0.05 < > 0.1" ="#a3cef1", "> 0.1"="white")) + 
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  geom_text(data = subset(df_heatmap_pigment_pv_long, Categorie != "> 0.1"), 
            aes(label = Arrow, color = Arrow_Color), size = 6) +
  scale_color_identity()+
  theme_test()+
  theme(axis.text.x = element_text(color = "black", face = "bold",  size = 18),
        axis.text.y = element_text(color = "black", face = "bold", size = 16),
        legend.position = "bottom")+
  annotate("segment", 
           x = 0.51, xend = 3.48,  
           y = 28, yend = 28,  
           colour = "black", 
           size = 1) +
  annotate("segment", 
           x = 3.51, xend = 6.5,  
           y = 28, yend = 28,  
           colour = "black", 
           size = 1) +
  annotate("text", 
           x = 2, y = 28.7, 
           label = "SB", 
           size = 8, 
           fontface = "bold", 
           colour = "black")+
  annotate("text", 
           x = 5, y = 28.7, 
           label = "WSB", 
           size = 8, 
           fontface = "bold", 
           colour = "black")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0.1)))

# FIGURE S11 ####

dosage_hp_df <- readxl::read_excel("/Users/alexandre/Desktop/Oxidative stress drive vertical migration MPB_metadata/SCRIPT/1_Data_formatting/data/Plasma_data/H2O2_calibration curve.xlsx")

# Standard calibration curve
cal_curve <- dosage_hp_df[is.na(dosage_hp_df$Time),]

linear_model <- lm(Absorbance ~ Concentration, data = cal_curve)

# Retrieve data
coeff <- coef(linear_model)
r2 <- summary(linear_model)$r.squared
equation <- paste0("y = ", format(coeff[2], digits =1, nsmall = 1),
                   "x + ", format(coeff[1], digits = 3, nsmall = 3),
                   "\nR² = ", format(r2, digits = 4, nsmall = 4))

# Plot
ggplot(cal_curve, aes(x = Concentration, y = Absorbance)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "#d62828") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
           label = equation, parse = FALSE, size = 5, color = "black") +
  labs(title = "",
       x = "H2O2 concentration (µM)",
       y = "Absorbance (λ409nm)") +
  theme_test()+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 14),
        axis.text.y = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15))


# TABLE S2 ####

pig_tab <- read_excel(here("data", "Pigments_data", "Pigment_diversity.xlsx"))

# Edit table
gt_table <- pig_tab %>% gt() %>% 
tab_style(
  style = cell_text(align = "left"),
  locations = cells_body(columns = everything())
) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything()))

gtsave(gt_table, filename = "ts2.png")

# TABLE S3 ####

df_Pparam <- df_Pparameters[, c("Name", "Time", "Day",
                                "Alpha", "Fq'/Fm' VLL", "rETRm obs", "F0", 
                                "NPQ1000", "YNPQ1000", "YNOm obs")]

df_Pparameters_1 <- subset(df_Pparam, Day == 1) # Day 1
df_Pparameters_2 <- subset(df_Pparam, Day == 2) # Day 2
df_Pparameters_3 <- subset(df_Pparam, Day == 3) # Day 3

# % differences ΔAfter-Before J1
df_Pparameters_1 <- df_Pparameters_1 %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = (((Alpha[Time == "After"]*100) / Alpha[Time == "Before"])-100),
    `Fq'/Fm' VLL` = (((`Fq'/Fm' VLL`[Time == "After"]*100) / `Fq'/Fm' VLL`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `F0` = (((`F0`[Time == "After"]*100) / `F0`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

# Mean sd J1
df_Pparameters_1 <- df_Pparameters_1 %>%
  mutate(Names = substr(Name, 1, nchar(Name) - 2)) %>%
  group_by(Names) %>%
  summarise(across(where(is.numeric), 
~ {if (all(is.na(.x))) NA_character_ else sprintf("%.1f (%.1f)", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE))}))

# % differences ΔAfter-Before J2
df_Pparameters_2 <- df_Pparameters_2 %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = (((Alpha[Time == "After"]*100) / Alpha[Time == "Before"])-100),
    `Fq'/Fm' VLL` = (((`Fq'/Fm' VLL`[Time == "After"]*100) / `Fq'/Fm' VLL`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `F0` = (((`F0`[Time == "After"]*100) / `F0`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

# Mean sd J2
df_Pparameters_2 <- df_Pparameters_2 %>%
  mutate(Names = substr(Name, 1, nchar(Name) - 2)) %>%
  group_by(Names) %>%
  summarise (across(where(is.numeric), 
~ {if (all(is.na(.x))) NA_character_ else sprintf("%.1f (%.1f)", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE))}))

# % differences ΔAfter-Before J3
df_Pparameters_3 <- df_Pparameters_3 %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `Alpha (α)` = (((Alpha[Time == "After"]*100) / Alpha[Time == "Before"])-100),
    `Fq'/Fm' VLL` = (((`Fq'/Fm' VLL`[Time == "After"]*100) / `Fq'/Fm' VLL`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `F0` = (((`F0`[Time == "After"]*100) / `F0`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

# Mean sd J3
df_Pparameters_3 <- df_Pparameters_3 %>%
  mutate(Names = substr(Name, 1, nchar(Name) - 2)) %>%
  group_by(Names) %>%
  summarise(across(where(is.numeric), 
~ {if (all(is.na(.x))) NA_character_ else sprintf("%.1f (%.1f)", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE))}))

# Transpositions
df_Pparameters_1<-t(df_Pparameters_1)
colnames(df_Pparameters_1) <- as.character(df_Pparameters_1[1, ])
df_Pparameters_1 <- df_Pparameters_1[-1, ]  
df_Pparameters_1 <- df_Pparameters_1[, c("SB_VLL","SB_HL","WSB_VLL","WSB_HL")]

df_Pparameters_2<-t(df_Pparameters_2)
colnames(df_Pparameters_2) <- as.character(df_Pparameters_2[1, ])
df_Pparameters_2 <- df_Pparameters_2[-1, ]  
df_Pparameters_2 <- df_Pparameters_2[, c("SB_VLL","SB_PL","WSB_VLL","WSB_PL")]

df_Pparameters_3<-t(df_Pparameters_3)
colnames(df_Pparameters_3) <- as.character(df_Pparameters_3[1, ])
df_Pparameters_3 <- df_Pparameters_3[-1, ]  
df_Pparameters_3 <- df_Pparameters_3[, c("SB_VLL","SB_HP","WSB_VLL","WSB_HP")]

# Print tables
df_Pparameters_1 <- as.data.frame(df_Pparameters_1, stringsAsFactors = FALSE)
df_Pparameters_1 <- tibble::rownames_to_column(df_Pparameters_1, var = "HL photosynthetic parameters")

gt_table <- df_Pparameters_1 %>%gt()%>%
  cols_label(
    `HL photosynthetic parameters` = "(A) HL photosynthetic parameters",
    SB_VLL = "SB VLL",
    SB_HL = "SB HL",
    WSB_VLL = "WSB VLL",
    WSB_HL = "WSB HL"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_column_labels(everything()),
      cells_body(columns = "HL photosynthetic parameters")))%>% 
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())) %>%
  sub_missing(missing_text = "")

df_Pparameters_2 <- as.data.frame(df_Pparameters_2, stringsAsFactors = FALSE)
df_Pparameters_2 <- tibble::rownames_to_column(df_Pparameters_2, var = "PL photosynthetic parameters")

gtsave(gt_table, filename = "ts3_1.png")

gt_table <- df_Pparameters_2 %>%gt()%>%
  cols_label(
    `PL photosynthetic parameters` = "(B) PL photosynthetic parameters",
    SB_VLL = "SB VLL",
    SB_PL = "SB PL",
    WSB_VLL = "WSB VLL",
    WSB_PL = "WSB PL"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_column_labels(everything()),
      cells_body(columns = "PL photosynthetic parameters")))%>% 
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything()))%>%
  sub_missing(missing_text = "")

gtsave(gt_table, filename = "ts3_2.png")

df_Pparameters_3 <- as.data.frame(df_Pparameters_3, stringsAsFactors = FALSE)
df_Pparameters_3 <- tibble::rownames_to_column(df_Pparameters_3, var = "HP photosynthetic parameters")

gt_table <- df_Pparameters_3 %>%gt()%>%
  cols_label(
    `HP photosynthetic parameters` = "(C) HP photosynthetic parameters",
    SB_VLL = "SB VLL",
    SB_HP = "SB HP",
    WSB_VLL = "WSB VLL",
    WSB_HP = "WSB HP"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_column_labels(everything()),
      cells_body(columns = "HP photosynthetic parameters")))%>% 
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())) %>%
  sub_missing(missing_text = "")

gtsave(gt_table, filename = "ts3_3.png")


# TABLE S4 ####
gdf_pigment_concentration <-gdf_pigment

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cDd = (Dd / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cDt = (Dt / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cV = (V / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cA = (A / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cZ = (Z / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration[, c("Treatment", "cDd", "cDt", "cV", "cA", "cZ")]

# Treatment name
gdf_pigment_concentration$Treatment <- gsub("_", " ", gdf_pigment_concentration$Treatment)

# Df structure
treatm_order <- c("VLL SB", "HL SB", "HP SB", "PL SB", "VLL WSB", "HL WSB", "HP WSB", "PL WSB")

# Mean sd calculation
tableS4 <- gdf_pigment_concentration %>%
  mutate(Treatment = factor(Treatment, levels = treatm_order)) %>%
  group_by(Treatment) %>%
  summarise(
    Diadinoxanthin = sprintf("%.1e (%.1e)", mean(cDd), sd(cDd)),
    Diatoxanthin = sprintf("%.1e (%.1e)", mean(cDt), sd(cDt)),
    Violaxanthin = sprintf("%.1e (%.1e)", mean(cV), sd(cV)),
    Antheraxanthin = sprintf("%.1e (%.1e)", mean(cA), sd(cA)),
    Zeaxanthin = sprintf("%.1e (%.1e)", mean(cZ), sd(cZ)))

# Edit table
gt_table <- tableS4 %>%
  gt() %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything()))

gtsave(gt_table, filename = "ts4.png")

# Table S5 ####
dosage_hp_df <- here("data", "Plasma_data", "H2O2_calibration curve.xlsx")
dosage_hp_df <- read_excel(dosage_hp_df)

dosage_hp <- dosage_hp_df[!is.na(dosage_hp_df$Time),]

# Subset
dosage_time_10 <- dosage_hp%>% filter(Time == 10)
dosage_time_5 <- dosage_hp%>% filter(Time == 5)

# Mean Sd Function
Add_stats <- function(df) {
  mean <- df %>% summarise(across(where(is.numeric), mean))
  sd <- df %>% summarise(across(where(is.numeric), sd))
  mean <- mean %>% mutate(Time = unique(df$Time), Note = "Mean")
  sd <- sd %>% mutate(Time = unique(df$Time), Note = "SD")
  df <- df %>% mutate(Note = "")
  bind_rows(df, mean, sd)
}

# Table 10 minutes
gt_table <- Add_stats(dosage_time_10) %>%
  gt() %>%
  tab_header(title = "PAW samples absorption after 10 minutes of plasma treatment") %>%
  cols_label(Note = "", Absorbance = "Absorbance (λ409nm)", Concentration = "H2O2 concentration (µM)") %>%
  cols_hide(columns = vars(Time)) %>%
  fmt_number(
    columns = vars(Absorbance),
    decimals = 3) %>%
  fmt_number(
    columns = vars(Concentration),
    decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = Note %in% c("Mean", "SD")))%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))


gtsave(gt_table, filename = "ts5_10.png")

# Table 5 minutes
gt_table <- Add_stats(dosage_time_5) %>%
  gt() %>%
  tab_header(title = "PAW samples absorption after 5 minutes of plasma treatment") %>%
  cols_label(Note = "", Absorbance = "Absorbance (λ409nm)", Concentration = "H2O2 concentration (µM)") %>%
  cols_hide(columns = vars(Time)) %>%
  fmt_number(
    columns = vars(Absorbance),
    decimals = 3) %>%
  fmt_number(
    columns = vars(Concentration),
    decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = Note %in% c("Mean", "SD")))%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))

gtsave(gt_table, filename = "ts5_5.png")

# TABLE S6 ####

pig_cal_curve <- readxl::read_excel("/Users/alexandre/Desktop/Oxidative stress drive vertical migration MPB_metadata/SCRIPT/1_Data_formatting/data/Pigments_data/DHI standard pigments/Calibration_curves.xlsx", sheet = 10)

# Edit table
gt_table <- pig_cal_curve %>%
  gt() %>%
  fmt_missing(
    columns = everything(),
    missing_text = ""
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything()))


gtsave(gt_table, filename = "ts6.png")