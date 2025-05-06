library(sf)

setwd(dirpath)
source("./R/builder.R")

# Use all catchments as seeds with a single area target
catchments_sample <- readRDS("data/catchments_sample.rds")

nghbrs <- neighbours(catchments_sample)
seed <- seeds(catchments_sf = catchments_sample, areatarget_value = 1000000000)
seed <- seed[1:10,]
builder(catchments_sf = catchments_sample, seeds = seed, neighbours = nghbrs, out_dir = file.path(dirpath, "Builder_output"))