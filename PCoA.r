# ==========================================
# 1. LOAD PACKAGES
# ==========================================
library(tidyverse)
library(reshape2)
library(vegan)

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
# 6. CALCULATE PCoA (AITCHISON)
# ==========================================
kraken_frac_t = t(kraken_frac)

# Aitchison distance 
dist_matrix = vegdist(kraken_frac_t, method = "aitchison")
pcoa_res = cmdscale(dist_matrix, k = 2, eig = TRUE)

# Extract percentage of explained variance
var_exp = round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)

pcoa_coords = as.data.frame(pcoa_res$points)
colnames(pcoa_coords) = c("PCoA1", "PCoA2")
pcoa_coords$SampleID = rownames(pcoa_coords)
pcoa_coords$SampleID = trimws(gsub("\\.", "-", pcoa_coords$SampleID))

# Merge coordinates with metadata
pcoa_meta = merge(pcoa_coords, meta, by = "SampleID")

# ==========================================
# 7. CALCULATE CENTROIDS FOR 4 GROUPS
# ==========================================
centroids_group = pcoa_meta %>%
  group_by(Gen, Group) %>%
  summarise(c_PCoA1 = mean(PCoA1), c_PCoA2 = mean(PCoA2), .groups = "drop")

# ==========================================
# 8. GENERATE COMBINED PCoA PLOT
# ==========================================
colors_group = c("FF" = "#CC79A7", 
                 "MF" = "#0072B2", 
                 "FS" = "#D55E00", 
                 "MS" = "#009E73")

pcoa_data_merged = merge(pcoa_meta, centroids_group, by = c("Gen", "Group"))

p_pcoa_combined = ggplot(pcoa_data_merged, aes(x = PCoA1, y = PCoA2, color = Group)) +
  stat_ellipse(level = 0.95, size = 0.6, show.legend = FALSE, alpha = 0.7) +
  geom_segment(aes(xend = c_PCoA1, yend = c_PCoA2), alpha = 0.5, size = 0.3, show.legend = FALSE) +
  geom_point(aes(x = c_PCoA1, y = c_PCoA2, fill = Group), 
             color = "white", shape = 21, size = 3.5, stroke = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", stroke = 0.8, size = 1.5) +
  scale_color_manual(values = colors_group) +
  scale_fill_manual(values = colors_group) +
  facet_wrap(~ Gen, ncol = 3) + 
  theme_bw() +
  labs(
    title = "PCoA of Microbiome by Phenotype Group\nAcross Generations",
    x = paste0("PCoA1 [", var_exp[1], "%]"),
    y = paste0("PCoA2 [", var_exp[2], "%]"),
    caption = "Groups: FF = Female Fighter, MF = Male Fighter, FS = Female Scrambler, MS = Male Scrambler"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 8.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = "black"), 
    strip.text = element_text(face = "bold", size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

print(p_pcoa_combined)

# Save plot to file
ggsave("/Users/janorzechowski/Library/Mobile Documents/com~apple~CloudDocs/Projekty - Magisterka/Praca Dyplomowa/Wstepne_Wyniki/PCoA_Aitchison_HighRes.png", 
       plot = p_pcoa_combined, 
       width = 11, height = 4.5, dpi = 300)

# ==========================================
# 9. GLOBAL PERMANOVA TEST (AITCHISON)
# ==========================================
# 1. Remove dots in sample names in Kraken matrix
rownames(kraken_frac_t) = trimws(gsub("\\.", "-", rownames(kraken_frac_t)))

# 2. Extract from 'meta' only those samples that are in Kraken 
# and arrange them in the same order (removes empty NA rows)
meta_aligned = meta[match(rownames(kraken_frac_t), meta$SampleID), ]

# 3. Calculate distance on Kraken
dist_all = vegdist(kraken_frac_t, method = "aitchison")

# 4. PERMANOVA using aligned metadata
perm_global = adonis2(dist_all ~ Sex * Type + Gen, data = meta_aligned, permutations = 999, by = "terms")
print(perm_global)

# ==========================================
# 10. GENERATE PCoA RESULTS TABLE (PERMANOVA)
# ==========================================
num_effects = nrow(perm_global) - 2

efekty = rownames(perm_global)[1:num_effects]
r2_surowe = perm_global$R2[1:num_effects]
p_surowe = perm_global[["Pr(>F)"]][1:num_effects]

beta_results_table = data.frame(
  Effect = efekty,
  R_squared_Percentage = round(r2_surowe * 100, 2),
  P_value = round(p_surowe, 4)
)

print(beta_results_table)