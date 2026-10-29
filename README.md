# Rhizoglyphus-Robini-Microbiome-Analysis
The data was obtained from a previous experiment conducted by the academic supervisor. The experiment involved the experimental evolution of Rhizoglyphus robini populations, aiming to investigate the impact of sexual selection on genome-wide genetic diversity. Whole-genome resequencing was carried out by the National Genomics Infrastructure (Uppsala, Sweden) using the Illumina Nova-Seq 6000 platform with an S4 flow cell to produce 2 × 150 bp reads (average 160.7 × 10<sup>6</sup>; range 130.7 × 10<sup>6</sup> - 189.9 × 10<sup>6</sup>).

## Repository Structure
* **`data/`**: Contains the sample metadata (`Indeks_Mite.csv`) and the final, preprocessed taxonomic count matrix (`kraken_family_matrix.csv`).
* **`scripts_R/`**: Contains R scripts for downstream statistical analysis and visualization.
  * `Alpha_Diversity.R`: Calculates Shannon and Simpson indices.
  * `PCoA.R`: Performs Aitchison distance transformation and PERMANOVA tests.
  * `Dynamic_Fractions.R`: Tracks family-level relative abundances across generations using Linear Mixed Models (LMM).
* **`pipeline.ipynb`**: A Jupyter Notebook detailing the upstream bioinformatics workflow used to generate the Kraken count matrix from raw FASTQ files.

## Downstream Analysis (R Modifications)
Several critical data transformations were implemented in the R downstream analysis to ensure statistical robustness:
* **Hominidae Removal:** Human DNA contamination (Hominidae family) was filtered out before calculating relative abundances to prevent artifactual skewing.
* **Statistical Modeling:** Linear Mixed Models (`~ Sex + Type + Gen + (1 | Line)`) were utilized for both Alpha Diversity and family dynamics.