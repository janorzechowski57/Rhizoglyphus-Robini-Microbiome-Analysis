# ==========================================
# 1. LOAD PACKAGES
# ==========================================
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(lme4)
library(lmerTest)

# ==========================================
# 2. LOAD FILES
# ==========================================
meta = read.csv("/Users/janorzechowski/Downloads/Magisterka_Analiza/Indeks_Mite.csv", sep = ";", stringsAsFactors = FALSE)
if(ncol(meta) == 1) {
  meta = read.csv("/Users/janorzechowski/Downloads/Magisterka_Analiza/Indeks_Mite.csv", sep = ",", stringsAsFactors = FALSE)
}
kraken = read.csv("/Users/janorzechowski/Downloads/Magisterka_Analiza/kraken_family_matrix.csv", stringsAsFactors = FALSE)

# ==========================================
# 3. CLEAN METADATA AND CREATE GROUPS
# ==========================================
colnames(meta)[1] = "SampleID"
meta$SampleID = trimws(meta$SampleID)
meta$SampleID = gsub("\\.", "-", meta$SampleID)

meta$Type = factor(meta$Type, levels = c("Fighter", "Scrambler"))
meta$Sex = factor(meta$Sex, levels = c("Female", "Male"))

# Create a numeric variable for the area plot and a categorical one for the bar plot
meta$Gen_num = as.numeric(gsub("[^0-9]", "", meta$Gen))
meta$Gen = factor(meta$Gen, levels = c("F1", "F12", "F29"))

# Combine variables and replace with abbreviations (FF, MF, FS, MS)
meta$Group = paste(meta$Sex, meta$Type, sep = " ")
meta$Group = factor(meta$Group, 
                    levels = c("Female Fighter", "Male Fighter", "Female Scrambler", "Male Scrambler"),
                    labels = c("FF", "MF", "FS", "MS"))

# ==========================================
# 4. PREPARE AND CLEAN KRAKEN DATA
# ==========================================
kraken_data = kraken[, -1]
rownames(kraken_data) = kraken[, 1]

# KEY CLEANING: Remove human DNA (Hominidae) before analysis
kraken_data = kraken_data[rownames(kraken_data) != "Hominidae", ]

# ==========================================
# 5. CALCULATE FRACTIONS xi = (ni + 1) / (N + m)
# ==========================================
m = nrow(kraken_data)
kraken_frac = apply(kraken_data, 2, function(n_i) {
  N = sum(n_i)
  x_i = (n_i + 1) / (N + m)
  return(x_i)
})

# ==========================================
# 6. SELECT TOP 7 FAMILIES
# ==========================================
mean_frac = rowMeans(kraken_frac)
top7_families = names(sort(mean_frac, decreasing = TRUE)[1:7])

final_data = kraken_frac[top7_families, ]
others_sum = 1 - colSums(final_data)
final_data = rbind(final_data, Others = others_sum)

# ==========================================
# 7. FORMATTING (MELT)
# ==========================================
final_data_melt = melt(final_data)
colnames(final_data_melt) = c("Family", "SampleID", "Fraction")
final_data_melt$SampleID = trimws(gsub("\\.", "-", final_data_melt$SampleID))

plot_data = merge(final_data_melt, meta, by = "SampleID")

# ==========================================
# 8. AVERAGE DATA (FOR BOTH PLOTS SIMULTANEOUSLY)
# ==========================================
# Group by both Gen and Gen_num to support both plot types
plot_data_group = plot_data %>%
  group_by(Group, Gen, Gen_num, Family) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

plot_data_group$Family = factor(plot_data_group$Family, levels = c(top7_families, "Others"))

# Define a common color palette for both plots
my_colors = c(brewer.pal(7, "Set2"), "grey80")

# ==========================================
# 9A. GENERATE AREA PLOT
# ==========================================
p_frac_area = ggplot(plot_data_group, aes(x = Gen_num, y = Fraction, fill = Family)) +
  geom_area(alpha = 0.85, color = "white", size = 0.2) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Group, ncol = 4) +
  theme_bw() +
  labs(
    title = "Microbiome Evolutionary Trends\nAcross Generations",
    x = "Generation",
    y = "Relative Abundance",
    caption = "Groups: FF = Female Fighter, MF = Male Fighter, FS = Female Scrambler, MS = Male Scrambler"
  ) +
  scale_x_continuous(breaks = c(1, 12, 29)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 8.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  )

# ==========================================
# 9B. GENERATE BAR PLOT
# ==========================================
p_frac_combined = ggplot(plot_data_group, aes(x = Group, y = Fraction, fill = Family)) +
  geom_col(width = 0.8, color = "white", size = 0.2) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Gen, ncol = 3) +
  theme_bw() +
  labs(
    title = "Microbiome Composition Dynamics\nAcross Generations", 
    x = "Phenotype Group",
    y = "Relative Abundance",
    caption = "Groups: FF = Female Fighter, MF = Male Fighter, FS = Female Scrambler, MS = Male Scrambler" 
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 8.5, color = "grey30"), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"), 
    axis.text.y = element_text(size = 11),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9), 
    legend.key.size = unit(0.4, "cm") 
  )

print(p_frac_area)
print(p_frac_combined)

# ==========================================
# 10. STATISTICAL ANALYSIS (LMM FOR FRACTIONS - ADDITIVE MODEL)
# ==========================================
# Logit transformation for proportions (0 < x < 1)
plot_data$Logit_Fraction = log(plot_data$Fraction / (1 - plot_data$Fraction))

# Loop to calculate the final LMM model for each bacterial family separately
lmm_results = plot_data %>%
  group_by(Family) %>%
  do({
    # Fit the mixed model with PLUS (without the non-significant interaction)
    mod = lmer(Logit_Fraction ~ Sex + Type + Gen + (1 | Line), data = .)
    
    # Extract ANOVA table (Type III) with p-values
    anova_tab = as.data.frame(anova(mod))
    anova_tab$Effect = rownames(anova_tab)
    
    # Return clean data frame
    anova_tab[, c("Effect", "F value", "Pr(>F)")]
  }) %>%
  ungroup()

# Rename columns for better readability
colnames(lmm_results) = c("Family", "Effect", "F_value", "P_value")

# Apply FDR correction (Benjamini-Hochberg) within each tested effect
lmm_results = lmm_results %>%
  group_by(Effect) %>%
  mutate(
    FDR_Adjusted_P = p.adjust(P_value, method = "BH"),
    Significant = ifelse(FDR_Adjusted_P < 0.05, "YES (*)", "NO")
  ) %>%
  ungroup() %>%
  arrange(Effect, FDR_Adjusted_P)

print(as.data.frame(lmm_results))