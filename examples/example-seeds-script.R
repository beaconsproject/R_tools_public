library(dplyr)
library(sf)
library(utils)

# Use all catchments as seeds with a single area target
setwd(dirpath)
source("./R/builder.R")

catchments_sample <- readRDS("data/catchments_sample.rds")
seeds(catchments_sf = catchments_sample, areatarget_value = 1000000000)

# Use column-based area targets
catchments_sample$area_target <- 1000000000
seeds(catchments_sf = catchments_sample, areatarget_col = "area_target")

# Filter based on a column value

catchments_sample %>%
  filter(intact == 1) %>%
  seeds(catchments_sf = ., areatarget_value = 1000000000)

# Filter using a spatial polygon
ref_poly <- data.frame(
  lon = c(-138.4, -138.1, -138.1, -138.4, -138.1, -138.1, -138, -138),
  lat = c(64.3, 64.3, 64.1, 64.1, 64.3, 64.1, 64.1, 64.3),
  Areatarget = c(rep(1000000000, 4), rep(2000000000, 4))
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  group_by(Areatarget) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  st_transform(st_crs(catchments_sample))

seeds(
  catchments_sf = catchments_sample,
  filter_polygon = ref_poly,
  areatarget_polygon = ref_poly,
  areatarget_polygon_col = "Areatarget"
)

# save seeds as csv
seed <- seeds(catchments_sf = catchments_sample, areatarget_value = 1000000000)
write.csv(seed, file=file.path(dirpath,"Builder_input/seeds.csv"), row.names=FALSE)
