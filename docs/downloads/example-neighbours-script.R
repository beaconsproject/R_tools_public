library(dplyr)
library(sf)
library(utils)

# Use all catchments as seeds with a single area target
setwd(dirpath)
source("./R/builder.R")

catchments_sample <- readRDS("data/catchments_sample.rds")
neighbours(catchments_sf = catchments_sample)

# save seeds as csv
neighbour <- neighbours(catchments_sf = catchments_sample)
write.csv(neighbour, file=file.path(dirpath,"Builder_input/neighbours.csv"), row.names=FALSE)
