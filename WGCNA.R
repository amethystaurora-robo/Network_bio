'''
This file takes the processed RNA counts and runs WGCNA to get modules.
It filters the modules based on the DEGs returned network_processing.ipynb.
It results in a csv with module colors ready to be visualized in Cytoscape.

'''

library(WGCNA)
library(grDevices)
library(RColorBrewer)

#load read counts and DEGs
read_counts_vst <- read.csv('rna_vst_proc.csv',header=FALSE)
degs <- read.csv('gestalt_dyngenie.csv')

#process data expression for WGCNA
datExpr <- read_counts_vst[-c(1, 2), ]  # Remove the first two rows
colnames(datExpr) <- datExpr[1, ]  # Set the first row as column names
datExpr <- datExpr[-1, ]            # Remove the first row
datExpr <- datExpr[, -1]             # Remove the first column (KEGG IDs with NANs)
head(datExpr)
rownames(datExpr) <- datExpr[, 1]   # Set the second column as row names
datExpr <- as.data.frame(t(datExpr))

#convert to matrix
datExpr <- as.matrix(datExpr)
head(datExpr)

#Pick soft threshold
power = pickSoftThreshold(datExpr, powerVector = seq(1, 20, by = 1), verbose = 5)

# Plotting Scale-Free topology fit index as a function of soft-thresholding power
par(mfrow = c(1, 2)) # Set up the plotting area

# Scale-Free Topology Fit
plot(power$fitIndices[, 1], power$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, R^2",
     type = "n", main = "Scale Free Topology Fit")

# Adding lines and points
text(power$fitIndices[, 1], power$fitIndices[, 2],
     labels = power$fitIndices[, 1], cex = 0.7, col = "red")

# Mean Connectivity Plot
plot(power$fitIndices[, 1], power$fitIndices[, 3],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")

# Adding lines and points
text(power$fitIndices[, 1], power$fitIndices[, 3],
     labels = power$fitIndices[, 1], cex = 0.7, col = "blue")

#Choose best soft threshold based on graphs
chosen_power <- 8

# Construct the adjacency matrix
adjacency = adjacency(datExpr, power = chosen_power)

# Create the topological overlap matrix (TOM)
TOM = TOMsimilarity(adjacency)

#get DEGs list from DynGENIE
#degs_list <- degs$KO

# Subset datExpr to only include genes in degs_list
#clean column names so they match correctly
#colnames(datExpr) <- gsub(";.*", "", colnames(datExpr))
#degs_subset = which(colnames(datExpr) %in% degs_list)
#datExpr_degs = datExpr[, degs_subset]
# Subset adjacency and TOM to include only the genes in degs_list
#adjacency_degs = adjacency[degs_subset, degs_subset]
#TOM_degs = TOM[degs_subset, degs_subset]

#Calculate the dissimilarity TOM
dissTOM = 1 - TOM

#Create the gene tree
geneTree_degs = hclust(as.dist(dissTOM), method = "average")

#Create the dynamicMods for the subset of genes
dynamicMods_degs = cutreeDynamic(dendro = geneTree_degs, distM = dissTOM, method = "tree", minClusterSize = 50)

moduleColors <- dynamicMods_degs  #numeric module assignments
uniqueModules <- unique(moduleColors)

# Combine colors from Set1, Set2, and Set3
palette1 <- brewer.pal(9, "Set1")    # Set1 has 9 colors
palette2 <- brewer.pal(8, "Set2")    # Set2 has 8 colors
palette3 <- brewer.pal(12, "Set3")   # Set3 has 12 colors

# Combine them into one palette
combinedPalette <- c(palette1, palette2, palette3)

# Generate additional colors by interpolation if needed
combinedPalette <- colorRampPalette(combinedPalette)(35)

# Assign these colors to modules
moduleColorAssignments <- combinedPalette[match(moduleColors, uniqueModules)]

#Plot the dendrogram and module colors for the subset of genes
plotDendroAndColors(geneTree_degs, moduleColorAssignments, "Module colors", dendroLabels = FALSE, hang = 0.03)

#Print out the unique module colors and corresponding modules
print(unique(moduleColorAssignments))

# 2. Create a data frame with gene names and their corresponding module colors
results_df = data.frame(Gene = geneNames, ModuleColor = dynamicMods_degs)
str(datExpr)
datExpr_cleaned <- subset(datExpr, select = -c(1, 2))  # Remove the first two columns for eigengenes
str(datExpr_cleaned)
length(datExpr_cleaned)
head(datExpr_cleaned)
# Assuming datExpr is your expression data and moduleColors are the colors assigned to each gene
MEs <- moduleEigengenes(datExpr_cleaned, colors = moduleColors)$eigengenes

# View the module eigengenes
head(MEs)

# Create a data frame mapping
moduleColorMapping <- data.frame(Module = uniqueModules, Color = colorPalette[uniqueModules + 1])  # +1 for 1-based indexing in R
str(datExpr)
# Check and Export the results to a CSV file
head(moduleColorMapping,50)
head(results_df,50)

final_df <- merge(results_df, moduleColorMapping, by.x = "ModuleColor", by.y = "Module", all.x = TRUE)
colnames(final_df)[colnames(final_df) == "Color"] <- "ModuleColorName"
head(final_df)
write.csv(final_df, file = "gene_module_colors.csv", row.names = FALSE)
final_df <- read.csv('gene_module_colors.csv')
print(length(final_df))
