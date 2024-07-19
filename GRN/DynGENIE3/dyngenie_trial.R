"""
This file runs the algorithm DynGENIE3 to infer a gene regulatory network.

It uses the wrapper doc dynGENIE3.R available from https://github.com/vahuynh/dynGENIE3
Input for the algorithm is pre-processed transcriptomic data for 3 conditions and a list of known 
transcription factors (TFs) for these genes

Input data is re-structured into lists, then the algorithm is run with the timepoints and gene 
expression lists and the list of known TFs as parameters

Output is a weighted adjacency matrix and a list of mRNA decay rates
These are written into csv files
"""

source('dynGENIE3.R')

library(reshape2)
library(igraph)
library(dplyr)

control <- read.expr.matrix('control_df.txt',form='rows.are.samples')
low <- read.expr.matrix('low_df.txt',form='rows.are.samples')
high <- read.expr.matrix('high_df.txt',form='rows.are.samples')
tf_terms <- read.csv('regulatory_genes.csv')
gene_order <- read.csv('gene_order.csv')

#run DynGENIE3 algorithm
run_dyngenie = function(df) {
  #re-structure data file into list of timepoints and expression levels
  time.points <- list(df[1,])
  TS.data <- list(df[2:nrow(df),])

  #run DynGENIE3
  res <- dynGENIE3(TS.data,time.points,regulators=c(tf_terms$regulatory),tree.method="ET",ntrees=50,K=5)

  #get weighted adjacency matrix
  link.list <- get.link.list(res$weight.matrix)

  # Extract the alphas values
  alphas <- res$alphas

  # Create a data frame with gene names and their corresponding alphas
  alphas_with_genes <- data.frame(gene = gene_order$Genes, alpha = alphas)

  #return weighted adjacency matrix and decay rates
  return (list(alphas_with_genes,link.list))
}

#run DynGENIE3 for 3 conditions
result_control <- run_dyngenie(control)
alphas_control <- result_control[[1]]
link.list_control <- result_control[[2]]

result_low <- run_dyngenie(low)
alphas_low <- result_low[[1]]
link.list_low <- result_low[[2]]

result_high <- run_dyngenie(high)
alphas_high <- result_high[[1]]
link.list_high <- result_high[[2]]

#write files of weighted adjacency matrix and decay rates
write.csv(link.list_control, "link_list_control.csv", row.names = FALSE)
write.csv(alphas_with_genes, 'alphas_control.csv', row.names=FALSE)

#write low files
write.csv(link.list_low, "link_list_low.csv", row.names = FALSE)
write.csv(alphas_with_genes_low, "alphas_low.csv", row.names = FALSE)

#write high files
write.csv(link.list_high, "link_list_high.csv", row.names = FALSE)
write.csv(alphas_with_genes_high, "alphas_high.csv", row.names = FALSE)

