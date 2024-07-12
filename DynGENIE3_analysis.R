control_weights <- read.csv('link_list_control.csv')
low_weights <- read.csv('link_list_low.csv')
high_weights <- read.csv('link_list_high.csv')
gene_order <- read.csv('gene_order.csv')
control_alphas <- read.csv('alphas_control.csv')
low_alphas <- read.csv('alphas_low.csv')
high_alphas <- read.csv('alphas_high.csv')

#get highest timepoint changes for a condition
get_top_alphas <- function(weights,alphas) {
  
  #merge alphas with gene order
  alphas_genes <- cbind(alphas,gene_order)
  
  #order by alphas, select top alphas
  control_ordered <- alphas_genes[order(alphas_genes$alpha, decreasing = TRUE),]
  top_5000_alphas <- control_ordered[control_ordered$alpha > 0.08,]
  
  # Rename the column "Gene" to "regulatory.gene"
  colnames(top_5000_alphas)[colnames(top_5000_alphas) == "Genes"] <- "target.gene"
  
  #merge alphas with main file
  
  weights_alphas <- merge(top_5000_alphas, weights, by = "target.gene")
  
  # Calculate node out degrees, sum weight for each transcription factor
  regulatory_sum <- aggregate(weight ~ regulatory.gene, data = weights_alphas, sum)
  colnames(regulatory_sum) <- c("regulatory.gene", "out degree")
  
  # Calculate node in degrees, sum weight for other genes
  target_sum <- aggregate(weight ~ target.gene, data = weights_alphas, sum)
  colnames(target_sum) <- c("target.gene", "in degree")
  
  # Merge degree for regulatory and target genes into alphas data frame
  degrees_target <- merge(weights_alphas, target_sum, by = "target.gene")
  degrees <- merge(degrees_target, regulatory_sum, by='regulatory.gene')
  
  return (degrees)
}


#get top alphas for control, low and high conditions
alphas_weights_control <- get_top_alphas(control_weights, control_alphas)
alphas_weights_low <- get_top_alphas(low_weights, low_alphas)
alphas_weights_high <- get_top_alphas(high_weights, high_alphas)

#write csv files for 3 conditions
write.csv(alphas_weights_control,'control_subset.csv')
write.csv(alphas_weights_low,'low_subset.csv')
write.csv(alphas_weights_high,'high_subset.csv')



