"""
This file finds genes which are differentially expressed under different conditions

It takes the input of merged alphas + decay weights for 3 conditions
Pairs are found between low and high dosed genes, all 3 conditions, only control, only low, and only high.
These filters are written to csv files. 

"""

library(dplyr)

control_8 <- read.csv('control_subset.csv')
low_8 <- read.csv('low_subset.csv')
high_8 <- read.csv('high_subset.csv')

#Get common genes which have alphas above 0.08 in all 3 conditions
# Extract target/regulatory gene pairs from each condition - control, low and high
# Create pairs_control with all columns from control_8
pairs_control <- control_8 %>%
  select(target.gene, regulatory.gene, everything()) %>%
  rename(target = target.gene, regulatory = regulatory.gene)

# Create pairs_low with all columns from low_8
pairs_low <- low_8 %>%
  select(target.gene, regulatory.gene, everything()) %>%
  rename(target = target.gene, regulatory = regulatory.gene)

# Create pairs_high with all columns from high_8
pairs_high <- high_8 %>%
  select(target.gene, regulatory.gene, everything()) %>%
  rename(target = target.gene, regulatory = regulatory.gene)

#Merge all 3 condition dfs to see genes which are common, result means these genes change as Daphnia grows/reproduces
#resulting file is control, x=low, y=high
common_pairs_cl <- merge(pairs_control, pairs_low, by=c('target', 'regulatory'))
common_pairs_clh <- merge(pairs_high, common_pairs_cl, by=c('target', 'regulatory'))

#merge pairs which are common in low and high doses but not control, this could be significant differences from ethoprophos
#resulting file is x=high, y=low
common_pairs_lh <- merge(pairs_high, pairs_low, by=c('target','regulatory'))
common_not_control <- anti_join(common_pairs_lh, pairs_control, by=c('target', 'regulatory'))

#find pairs only in control
pairs_only_control <- setdiff(pairs_control,pairs_low)
pairs_control_only <- setdiff(pairs_only_control,pairs_high)

#find pairs only in low
pairs_only_low <- setdiff(pairs_low,pairs_control)
pairs_low_only <- setdiff(pairs_only_low,pairs_high)

#find pairs only in high
pairs_only_high <- setdiff(pairs_high,pairs_low)
pairs_high_only <- setdiff(pairs_only_high,pairs_control)

write.csv(common_not_control,'dosed_genes_8.csv')
write.csv(common_pairs_clh,'growth_reproduction_genes_8.csv')
write.csv(pairs_control_only,'control_only_genes_8.csv')
write.csv(pairs_low_only,'low_only_genes_8.csv')
write.csv(pairs_high_only,'high_only_genes_8.csv')
