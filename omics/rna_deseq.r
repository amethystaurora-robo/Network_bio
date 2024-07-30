"""
This file runs DeSeq2 using Log2 fold change shrinkage.
It takes the metadata and raw data from rna_preproc.r.
This file will return a list of DEGs and their metrics.
It also generates PCA plots.

"""

# Loading the Libraries and Setting the Working Path
# ---------------------------------------------------
rm(list = ls())

library(DESeq2)
library(ashr)
library(RColorBrewer)


# Loading the Sample Sheet and Read Counts
# ---------------------------------------------------
# Load sample sheet
sample_sheet_file <- 'metadata_raw_filtered.csv'
sample_sheet <- read.table(sample_sheet_file, sep = ',', header = T, row.names = 1, stringsAsFactors = T, check.names = F, comment.char = '')

head(sample_sheet)
dim(sample_sheet)


# Load read counts
read_counts_file <- 'rna_raw_filtered.csv'
read_counts <- read.table(read_counts_file, sep = ',', header = T, row.names = 1, check.names = F, comment.char = '')

head(read_counts)
dim(read_counts)

read_counts <- as.matrix(read_counts)


# Check the order of sample sheet and read counts
stopifnot(row.names(sample_sheet) == colnames(read_counts))
stopifnot(!any(is.na(row.names(sample_sheet))))
stopifnot(!any(is.na(colnames(read_counts))))

dim(read_counts)

# Subset treatment group
# ---------------------------------------------------
sub_ids <- sample_sheet$treatment %in% c('low', 'high','control')
sub_sample_sheet <- sample_sheet[sub_ids,]
sub_read_counts <- read_counts[,sub_ids]

sub_sample_sheet  <- droplevels(sub_sample_sheet)
sub_sample_sheet$time_point <- factor(sub_sample_sheet$time_point, 
                                  levels = c('1H', '2H', '6H', '12H', '24H', '4D', '5D', '6D', '7D'))

sub_sample_sheet$time_point <- relevel(sub_sample_sheet$time_point, ref = '1H')

sub_sample_sheet$treatment <- relevel(sub_sample_sheet$treatment, ref = 'control') 

levels(sub_sample_sheet$time_point)

# DESeq2 on individual timepoint
# ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  sub_read_counts, colData = sub_sample_sheet,
  design = ~ time_point + treatment + time_point:treatment
) 

for (tp in levels(sub_sample_sheet$time_point)){
  sub_dds <- dds[dds$time_point == tp, ]
  sub_dds <- DESeq(sub_dds)
}




dds <- DESeqDataSetFromMatrix(
    raw_counts, colData = sample_sheet,
    # design = ~ treatment * timepoint
    design = ~ time_point + treatment + time_point:treatment
)

# Reorder time_point levels
sample_sheet$time_point <- factor(sub_sample_sheet$time_point, 
                                      levels = c('1H', '2H', '6H', '12H', '24H', '4D', '5D', '6D', '7D'))

# Check the new levels
levels(sample_sheet$time_point)
levels(sub_sample_sheet$time_point)

# DESeq2 fitting
# dds <- DESeq(dds)
ddsLRT = DESeq(dds, test='LRT', reduced= ~ time_point + treatment)
# plotDispEsts(dds)


resultsNames(ddsLRT)[-1]
resultsNames(ddsLRT)
for (tp in c('1H','2H','6H','12H','24H','4D','5D','6D','7D')){
  for (tr in c('low','high')){
    res.name <- paste0(tp,'_',tr,'vC')
    cat('running ', res.name, '\n')
    cat('------------------------------------\n')
    
    if (tp == '1H'){
      res.comp <- paste0('treatment_',tr,'_vs_control')
      res <- results(ddsLRT, name=res.comp, test='Wald')
    } else {
      res.contrast <- list(c(paste0('treatment_',tr,'_vs_control'), paste0('time_point',tp,'.treatment',tr)))
      res <- results(ddsLRT, contrast=res.contrast, test='Wald')
    }
      
    # number of DEGs (adjusted p-values < 0.05)
    cat(paste0(res.name, ' significant p-values = ', sum(res$pvalue < 0.05, na.rm = T), '\n'))
    cat(paste0(res.name, ' significant p-adjust = ', sum(res$padj < 0.05, na.rm = T), '\n'))

    # lfcs
    if (tp == '1H'){
      lfc_res <- lfcShrink(dds = ddsLRT, coef = res.comp, res = res)
    } else {
      lfc_res <- lfcShrink(dds = ddsLRT, contrast = res.contrast, res = res, type = 'ashr')
    }

    # export
    write.csv(lfc_res, paste0(res.name, '_deseq_res.csv'))
    
    cat('\n')
  }
}

# Obtain normalised read counts
dds <- estimateSizeFactors(dds) # if not fitted the full model
norm_counts <- counts(dds, normalized = T)
write.csv(norm_counts, file = 'rna_norm_counts.csv')

# Obtain vst transformed read counts
vds <- varianceStabilizingTransformation(dds, blind = TRUE) # if not fitted the full model
# vds <- varianceStabilizingTransformation(dds) #, fitType = 'local')
vst_counts <- as.matrix(assay(vds))
write.csv(vst_counts,  file = 'rna_vst_counts.csv')

# Obtain rlog transformed read counts
# rds <- rlog(dds)
# log_counts <- as.matrix(assay(rds))
# write.csv(log_counts[genes,],  file = paste0('rna_rlog.csv'))


# PCA
# ---------------------------------------------------
# PCA plot (pick the outliers and re-fit the DESeq model)
pca.plot <- function(read.counts, classes, 
                     comps = c(1, 2), ntop = min(500, nrow(read.counts)), standard = T,
                     col = c('lightblue', 'orange', 'MediumVioletRed', 'SpringGreen')){
  top_index <- order(rowVars(as.matrix(read.counts), useNames = FALSE), decreasing = TRUE)[1:ntop]
  pca <- prcomp(scale(t(read.counts[top_index,]), center = standard, scale = standard))
  prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
  
  pca_comps <- pca$x[,comps]
  prp_comps <- round(prp[comps], 2)

  df <- data.frame(pc1 = pca_comps[,1], pc2 = pca_comps[,2], condition = classes)
  p  <- ggplot(df, aes(x = pc1, y = pc2, color = condition)) + 
        geom_point(size = 3) + 
        labs(title = paste0('Principal Component Analysis - Axes ', comps[1] , ' and ', comps[2]), 
             x = paste0('PC', comps[1], ' (', prp_comps[1], '%)'), 
             y = paste0('PC', comps[2], ' (', prp_comps[2], '%)')) + 
        # geom_text(label = colnames(read.counts), vjust = 0, nudge_y = 2) +
        scale_color_manual(values = col)
  return(p)
}

# display.brewer.all()

# PCA plot on read counts
sids <- sample_sheet$treatment != 'control'

groups <- as.factor(sample_sheet$time_point)
groups <- as.factor(paste(sample_sheet$time_point, sample_sheet$treatment))
counts <- vst_counts
print(groups)

p <- pca.plot(counts, groups, comps = c(1,2), ntop = dim(counts)[1],
              col = rep(brewer.pal(min(max(3, nlevels(groups)),12), 'Paired'), 10))

p <- pca.plot(counts[,sids], groups[sids], comps = c(1,2), ntop = 2000,
              col = rep(brewer.pal(min(max(3, nlevels(groups)),12), 'Paired'), 10))
# p
ggplotly(p)


