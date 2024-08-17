"""
This file finds genes which are differentially expressed under different conditions

It takes the input of merged alphas + decay weights for 3 conditions
Genes under treatment conditions (low and high) are merged and annotated with the treatment.

"""

library(dplyr)

control_8 <- read.csv('control_subset.csv')
low_8 <- read.csv('low_subset.csv')
high_8 <- read.csv('high_subset.csv')

library(dplyr)

control_8 <- read.csv('control_subset.csv')
low_8 <- read.csv('low_subset.csv')
high_8 <- read.csv('high_subset.csv')

# Create pairs with consistent column names and add source column
pairs_low <- low_8 %>%
  select(target.gene, regulatory.gene, everything()) %>%
  rename(target = target.gene, regulatory = regulatory.gene) %>%
  mutate(source = 'low')

pairs_high <- high_8 %>%
  select(target.gene, regulatory.gene, everything()) %>%
  rename(target = target.gene, regulatory = regulatory.gene) %>%
  mutate(source = 'high')

# Combine low and high, adding source information
combined <- bind_rows(pairs_low, pairs_high)

# Exclude pairs that are present in control
pairs_control <- control_8 %>%
  select(target.gene, regulatory.gene, everything()) %>%
  rename(target = target.gene, regulatory = regulatory.gene)

result <- anti_join(combined, pairs_control, by=c('target', 'regulatory'))

head(result)
# Save the result to a CSV file
write.csv(result, 'genes_in_low_or_high_with_source.csv', row.names = FALSE)

