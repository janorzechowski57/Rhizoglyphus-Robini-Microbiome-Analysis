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
# 6. GENERATE 2 BOXPLOT + JITTER PLOTS
# ==========================================
# Define two colors for sex
colors_sex = c("Female" = "#CC79A7", "Male" = "#0072B2")

# -- 1. SHANNON --
p_shannon = ggplot(plot_data, aes(x = Type, y = Shannon, fill = Sex, shape = Sex)) +
  facet_wrap(~ Gen) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.45, position = position_dodge(width = 0.65)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.65), size = 2.5, color = "black", alpha = 0.7) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = colors_sex) +
  scale_shape_manual(values = c("Female" = 21, "Male" = 22)) +
  labs(
    title = "Shannon Diversity by Phenotype Group\nAcross Generations", 
    x = "Morph Treatment", 
    y = "Shannon Diversity Index",
    fill = "Sex",
    shape = "Sex"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.margin = margin(t = 0), 
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 5)),
    axis.text.x = element_text(size = 9.5, face = "bold", color = "black"), 
    axis.text.y = element_text(size = 11),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# -- 2. SIMPSON --
p_simpson = ggplot(plot_data, aes(x = Type, y = Simpson, fill = Sex, shape = Sex)) +
  facet_wrap(~ Gen) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.45, position = position_dodge(width = 0.65)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.65), size = 2.5, color = "black", alpha = 0.7) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = colors_sex) +
  scale_shape_manual(values = c("Female" = 21, "Male" = 22)) +
  labs(
    title = "Simpson Diversity by Phenotype Group\nAcross Generations", 
    x = "Morph Treatment", 
    y = "Simpson Diversity Index",
    fill = "Sex",
    shape = "Sex"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.margin = margin(t = 0), 
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 5)),
    axis.text.x = element_text(size = 9.5, face = "bold", color = "black"), 
    axis.text.y = element_text(size = 11),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
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