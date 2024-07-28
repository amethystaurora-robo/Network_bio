"""
This file is used for parameter tuning of the model DynGENIE3. 
The control file is used for parameter tuning. 
DynGENIE3 has been run on the control file with known transcription factors to provide a
ground truth to compare the model without known transcription factors. Regulatory 
relationships which have been implicated with provided known transcription factors are
compared to running the model without knowledge of these regulatory relationships.

In this file, 7 values of k (randomly chosen candidate regulators at each node of a tree) and 5
values of ntrees (number of trees per ensemble) are tested against the ground truth file.
Each combination is tried and tested against ground truth. Three metrics - precision, recall and f1 
are obtained. The parameters with the best f1 score are selected. In this case, a k of 5 and ntrees of
50 has the highest f1 score - 0.947.

"""
source('dynGENIE3.R')

library(dplyr)

control <- read.expr.matrix('control_df.txt',form='rows.are.samples')

#re-format timepoints and data for DynGENIE3
time.points_control <- list(control[1,])
TS.data_control <- list(control[2:nrow(control),])
tf_terms <- read.csv('regulatory_genes.csv')


res_control <- dynGENIE3(TS.data_control,time.points_control,regulators=c(tf_terms$regulatory))
#get weighted adjacency matrix
link.list_control <- get.link.list(res_control$weight.matrix)
num_top_quarter <- ceiling(0.25 * nrow(link.list_control))
top_quantile_known_TFs <- link.list_control[1:num_top_quarter, ]
write.csv(top_quantile_known_TFs,'top_quantile_known_TFS')
head(top_quantile_known_TFs)
#run Dyngenie3 with known TFs, get link list in order of weights, find length and extract top 1/4 (upper quantile)
#use list of upper quantile, loop and compare with dyngenie3 

param_grid <- expand.grid(
  K = c(3, 5, 10, 20, 50, 100, 150),        # Expanded range of candidate regulators
  n_trees = c(50, 100, 200, 300, 500)  # Number of trees
)

# Function to evaluate the network performance
evaluate_grn <- function(predicted_grn, ground_truth) {
  predicted_set <- paste(predicted_grn$regulatory.gene, predicted_grn$target.gene, sep = "->")
  ground_truth_set <- paste(ground_truth$regulatory.gene, ground_truth$target.gene, sep = "->")
  
  true_positives <- sum(predicted_set %in% ground_truth_set)
  false_positives <- sum(!predicted_set %in% ground_truth_set)
  false_negatives <- sum(!ground_truth_set %in% predicted_set)
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(f1_score)
}

# Initialize a results data frame
results <- data.frame(K = integer(), n_trees = integer(), f1_score = numeric())

# Perform grid search
for (i in 1:nrow(param_grid)) {
  K <- param_grid$K[i]
  n_trees <- param_grid$n_trees[i]
  
  # Run dynGENIE3 with the specified parameters
  res_unknown_TFs <- dynGENIE3(TS.data_control, time.points_control, K = K, ntrees = n_trees)
  link.list_unknown_TFs <- get.link.list(res_unknown_TFs$weight.matrix)
  top_quantile_unknown_TFs <- link.list_unknown_TFs[1:num_top_quarter, ]
  print(n_trees)
  
  # Evaluate the GRN
  if (!is.null(top_quantile_known_TFs)) {
    evaluation <- evaluate_grn(top_quantile_unknown_TFs,top_quantile_known_TFs)
    f1_score <- evaluation
  } else {
    # If ground truth is not available, just store NA for evaluation
    f1_score <- NA
  }
  
  # Store the results
  results <- rbind(results, data.frame(K = K, n_trees = n_trees, f1_score = f1_score))
}

# Remove the first empty row (initialized row)
results <- results[-1, ]

# Find the best parameters based on F1-score
best_params <- results[which.max(results$f1_score), ]
write.csv(results,'results.csv')
