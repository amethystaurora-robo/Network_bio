"""
This file uses pre-processed transcriptomic data and metadata from omics/deseq_pre-processing.ipynb
It uses Limma to obtain a list of DEGs between experimental conditions
and visualizes these genes with volcano and MA plots and shows the distribution of p-values for DEGs.
"""

library(limma)
library(edgeR) 
library(ggplot2)
library(pheatmap)

rna_counts <- read.csv("raw_counts_cleaned.csv", row.names = 1)
metadata_df <- read.csv("metadata_rna.csv", row.names = 1)

#convert metadata to numeric for limma 
metadata_df$conditions <- factor(metadata_df$conditions,levels = c('control', 'low', 'high'))
metadata_df$timepoints <- factor(metadata_df$timepoints, levels = c('1H', '2H', '6H', '12H', '24H', '4D', '5D', '6D', '7D'))

# Create design matrix
design <- model.matrix(~ 0 + conditions*timepoints, data = metadata_df)
colnames(design) <- make.names(colnames(design))

# Fit linear model
fit <- lmFit(rna_counts, design)

# Define contrasts for each comparison
contrasts <- makeContrasts(
  Control_vs_Low = conditionslow - conditionscontrol,
  Control_vs_High = conditionshigh - conditionscontrol,
  Low_vs_High = conditionshigh - conditionslow,
  levels = design
)

analyze_contrasts <- function(results_name,coef_string,volcano_title,MA_title,pvalue_title) {

  # Apply contrasts to the fitted model
  fit2 <- contrasts.fit(fit, contrasts)

  # Apply eBayes to get moderated t-statistics and p-values
  fit2 <- eBayes(fit2)

  # Extract results for specific contrasts
  results_name <- topTable(fit2, coef=coef_string)

  # Volcano plot for Control vs Low at 1H
  ggplot(results_name, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(aes(color=adj.P.Val < 0.05)) +
    labs(title=volcano_title,
         x="Log Fold Change",
         y="-Log10 Adjusted P-Value") +
    theme_minimal()

  # MA plresults_Control_vs_Low_24H# MA plot for Control vs Low at 1H
  plotMA(fit2, coef=coef_string, main=MA_title)

  #look at distribution of pvalues
  hist(results_name$P.Value, breaks=50, main=pvalue_title, xlab="P-value")
}

results_Control_vs_Low <- analyze_contrasts(
  coef_string = "Control_vs_Low",
  volcano_title = "Volcano Plot: Control vs Low",
  MA_title = "MA Plot: Control vs Low",
  pvalue_title = "P-value Distribution: Control vs Low"
)

results_Control_vs_High <- analyze_contrasts(
  coef_string = "Control_vs_High",
  volcano_title = "Volcano Plot: Control vs High",
  MA_title = "MA Plot: Control vs High",
  pvalue_title = "P-value Distribution: Control vs High"
)

results_Low_vs_High <- analyze_contrasts(
  coef_string = "Low_vs_High",
  volcano_title = "Volcano Plot: Low vs High",
  MA_title = "MA Plot: Low vs High",
  pvalue_title = "P-value Distribution: Low vs High"
)
