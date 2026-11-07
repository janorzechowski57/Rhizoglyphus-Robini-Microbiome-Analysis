# ==========================================
# 1. LOAD PACKAGES
# ==========================================
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)

# ==========================================
# 2. LOAD FILES
# ==========================================
meta = read.csv("/Users/janorzechowski/Documents/GitHub/Rhizoglyphus-Robini-Microbiome-Analysis/data/Indeks_Mite.csv", sep = ";", stringsAsFactors = FALSE)
if(ncol(meta) == 1) {
  meta = read.csv("/Users/janorzechowski/Documents/GitHub/Rhizoglyphus-Robini-Microbiome-Analysis/data/Indeks_Mite.csv", sep = ",", stringsAsFactors = FALSE)
}
kraken = read.csv("/Users/janorzechowski/Documents/GitHub/Rhizoglyphus-Robini-Microbiome-Analysis/data/kraken_family_matrix.csv", stringsAsFactors = FALSE)

# ==========================================
# 3. CLEAN METADATA AND CREATE GROUPS
# ==========================================
colnames(meta)[1] = "SampleID"
meta$SampleID = trimws(meta$SampleID)
meta$SampleID = gsub("\\.", "-", meta$SampleID)

meta$Type = factor(meta$Type, levels = c("Fighter", "Scrambler"))
meta$Sex = factor(meta$Sex, levels = c("Female", "Male"))
meta$Gen = factor(meta$Gen, levels = c("F1", "F12", "F29"))
meta$Line = factor(meta$Line)

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
# 5. CALCULATE ALPHA DIVERSITY INDICES (SHANNON AND SIMPSON)
# ==========================================
# Transpose cleaned Kraken matrix (samples as rows, families as columns)
kraken_t = as.data.frame(t(kraken_data))
kraken_t[] = lapply(kraken_t, as.numeric)
kraken_t$SampleID = trimws(gsub("\\.", "-", rownames(kraken_t)))

# Merge metadata with cleaned Kraken counts
merged_data = merge(meta, kraken_t, by = "SampleID")
meta_aligned = merged_data[, 1:ncol(meta)]
kraken_aligned = merged_data[, (ncol(meta) + 1):ncol(merged_data)]

# Calculate diversity indices
clean_shannon = diversity(kraken_aligned, index = "shannon")
clean_simpson = diversity(kraken_aligned, index = "simpson")

plot_data = data.frame(
  SampleID = meta_aligned$SampleID,
  Group = meta_aligned$Group, 
  Line = meta_aligned$Line,
  Sex = meta_aligned$Sex, 
  Type = meta_aligned$Type,
  Gen = meta_aligned$Gen,
  Shannon = clean_shannon,
  Simpson = clean_simpson
)

# ==========================================
# 6. GENERATE SUMMARY PLOTS (NATURE ECOLOGY STYLE)
# ==========================================
# Calculate means and standard deviations including Sex
plot_summary = plot_data %>%
  group_by(Type, Gen, Sex) %>%
  summarise(
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Simpson = mean(Simpson),
    SD_Simpson = sd(Simpson),
    .groups = 'drop'
  )

colors_type = c("Fighter" = "#24A385", "Scrambler" = "#DE7422")
shapes_sex = c("Female" = 21, "Male" = 22) 

# Position dodge to separate points slightly
pd = position_dodge(0.7)

# -- 1. SHANNON --
p_shannon = ggplot() +
  # Raw data (small points) - FIX: Added group = interaction(Type, Sex) to strictly align with summary points
  geom_jitter(data = plot_data, aes(x = Gen, y = Shannon, color = Type, shape = Sex, group = interaction(Type, Sex)), 
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), 
              size = 1.5, alpha = 0.4, show.legend = FALSE) +
  # Error bars
  geom_errorbar(data = plot_summary, aes(x = Gen, ymin = Mean_Shannon - SD_Shannon, 
                                         ymax = Mean_Shannon + SD_Shannon, color = Type, group = interaction(Type, Sex)), 
                position = pd, width = 0.15, size = 1, show.legend = TRUE) +
  # Mean points
  geom_point(data = plot_summary, aes(x = Gen, y = Mean_Shannon, color = Type, shape = Sex, group = interaction(Type, Sex)), 
             position = pd, size = 4, stroke = 1.5, fill = "white", show.legend = TRUE) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = colors_type, 
                     labels = c("Fighter" = "Fighter", "Scrambler" = "Scrambler"),
                     guide = guide_legend(override.aes = list(shape = 21, fill = "white", size = 4, stroke = 1.5, linetype = 1))) +
  scale_shape_manual(values = shapes_sex) +
  labs(
    title = "Shannon Diversity by Phenotype Group\nAcross Generations", 
    y = NULL, 
    x = "Generation",
    color = "Morph",
    shape = "Sex"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# -- 2. SIMPSON --
p_simpson = ggplot() +
  # Raw data (small points) - FIX: Added group = interaction(Type, Sex)
  geom_jitter(data = plot_data, aes(x = Gen, y = Simpson, color = Type, shape = Sex, group = interaction(Type, Sex)), 
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), 
              size = 1.5, alpha = 0.4, show.legend = FALSE) +
  # Error bars
  geom_errorbar(data = plot_summary, aes(x = Gen, ymin = Mean_Simpson - SD_Simpson, 
                                         ymax = Mean_Simpson + SD_Simpson, color = Type, group = interaction(Type, Sex)), 
                position = pd, width = 0.15, size = 1, show.legend = TRUE) +
  # Mean points
  geom_point(data = plot_summary, aes(x = Gen, y = Mean_Simpson, color = Type, shape = Sex, group = interaction(Type, Sex)), 
             position = pd, size = 4, stroke = 1.5, fill = "white", show.legend = TRUE) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = colors_type, 
                     labels = c("Fighter" = "Fighter", "Scrambler" = "Scrambler"),
                     guide = guide_legend(override.aes = list(shape = 21, fill = "white", size = 4, stroke = 1.5, linetype = 1))) +
  scale_shape_manual(values = shapes_sex) +
  labs(
    title = "Simpson Diversity by Phenotype Group\nAcross Generations", 
    y = NULL,
    x = "Generation",
    color = "Morph",
    shape = "Sex"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p_shannon)
print(p_simpson)

# ==========================================
# 7. STATISTICAL ANALYSIS (LMM FOR ALPHA DIVERSITY - ADDITIVE MODEL)
# ==========================================
# 7A. MODEL FOR SHANNON INDEX
model_shannon = lmer(Shannon ~ Sex + Type + Gen + (1 | Line), data = plot_data)

cat("\n--- LMM RESULTS FOR SHANNON INDEX ---\n")
anova_shannon = anova(model_shannon)
print(anova_shannon)

# 7B. MODEL FOR SIMPSON INDEX
model_simpson = lmer(Simpson ~ Sex + Type + Gen + (1 | Line), data = plot_data)

cat("\n--- LMM RESULTS FOR SIMPSON INDEX ---\n")
anova_simpson = anova(model_simpson)
print(anova_simpson)

# ==========================================
# 8. VERIFY MODEL ASSUMPTIONS (SHAPIRO-WILK TEST)
# ==========================================
residuals_shannon = resid(model_shannon)
residuals_simpson = resid(model_simpson)

cat("\nNormality test of model residuals (Shannon):\n")
print(shapiro.test(residuals_shannon))

cat("\nNormality test of model residuals (Simpson):\n")
print(shapiro.test(residuals_simpson))