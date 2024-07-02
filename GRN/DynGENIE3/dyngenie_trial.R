source('dynGENIE3.R')

control <- read.expr.matrix('control_df.txt',form='rows.are.samples')
low <- read.expr.matrix('low_df.txt',form='rows.are.samples')
high <- read.expr.matrix('high_df.txt',form='rows.are.samples')
head(control)
TS1 <- read.expr.matrix('time_series_1.txt',form='rows.are.samples')
TS2 <- read.expr.matrix('time_series_2.txt',form='rows.are.samples')
TS3 <- read.expr.matrix('time_series_3.txt',form='rows.are.samples')

subset <- read.expr.matrix('subset_df.txt',form='rows.are.samples')

#re-format example data
time.points <- list(control[1,],low[1,],high[1,]) 
TS.data <- list(control[2:nrow(control),],low[2:nrow(low),],high[2:nrow(high),])
#genes are in rows, timepoints in columns, each matrix represents a time series experiment
res <- dynGENIE3(TS.data,time.points)


#download, load example data
#run, print to see what exact format we need
#convert to correct format, run on my data
#clean up code, organize files, upload to git, upload to OneDrive

