## code to prepare `builder_catchments_sample` dataset goes here

library(sf)
library(dplyr)
library(devtools)

# load catchments from kba project (just using this dataset because its synced on my laptop, could also load full dataset from gisdata)
#catchments <- sf::st_read("C:/Users/MAEDW7/Dropbox (BEACONs)/RENR491 Capstone 2022/gisdata/catchments/YRW_catch50K.shp")
catchments <- sf::st_read("C:/Users/mehou10/BEACONs Dropbox/Mélina Houle/gisdata/catchments/boreal_vPB25/boreal_vPB25.shp")

catchments_sample <- catchments %>%
  #dplyr::filter(FDAHUC8 == "09EA") %>%
  dplyr::filter(FDA_M == "09EA") %>%
  #dplyr::select(CATCHNUM, SKELUID, STRAHLER, ORDER1, ORDER2, ORDER3, BASIN, Area_Land, Area_Water, Area_Total, STRMLEN_1, FDAHUC8, ZONE, MDA, Isolated, intact) %>%
  dplyr::select(CATCHNUM, SKELUID, STRAHLER, ORDER1, ORDER2, ORDER3, BASIN, Area_land, Area_water, Area_total, length_m, FDA_M, MDAzone, Isolated, CA2010) %>%
  sf::st_snap(x = ., y = ., tolerance = 0.1)

# change names to match previously used catchments
names(catchments_sample)[names(catchments_sample) == "Area_land"] <- "Area_Land"
names(catchments_sample)[names(catchments_sample) == "Area_water"] <- "Area_Water"
names(catchments_sample)[names(catchments_sample) == "Area_total"] <- "Area_Total"
names(catchments_sample)[names(catchments_sample) == "FDA_M"] <- "FDA"
names(catchments_sample)[names(catchments_sample) == "length_m"] <- "STRMLEN"
names(catchments_sample)[names(catchments_sample) == "MDAzone"] <- "MDA"
names(catchments_sample)[names(catchments_sample) == "CA2010"] <- "intact"

catchments_sample$CATCHNUM <- as.integer(catchments_sample$CATCHNUM)
catchments_sample$SKELUID <- as.integer(catchments_sample$SKELUID)
catchments_sample$ORDER2 <- as.numeric(catchments_sample$ORDER2)

catchments_sample$MDA <- as.character(catchments_sample$MDA)
catchments_sample$BASIN <- as.character(catchments_sample$BASIN)
catchments_sample$Isolated <- as.integer(catchments_sample$Isolated)

saveRDS(catchments_sample, file = "data/catchments_sample.rds")


# get existing reserves in FDA 09EA
temp <- file.path(tempdir(), "ProtectedConservedArea_2023.zip")
download.file("https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fprotectrestore%2Fcanadian-protected-conserved-areas-database%2FDatabases%2FProtectedConservedArea_2023.zip", temp) # unzip manually in temp file
reserves <- st_read(file.path(tempdir(), "ProtectedConservedArea_2023.gdb"), layer = "ProtectedConservedArea_2023")

reserves <- st_transform(reserves, st_crs(catchments_sample))
names(reserves)[names(reserves) == "Shape"] <- "geometry"
st_geometry(reserves)="geometry"

fda <- catchments_sample %>%
  summarise(geometry = st_union(geometry))

reserves <- reserves[reserves$NAME_E == "Tombstone",]

reserves_fda <- reserves %>%
  st_intersection(fda) %>%
  st_cast("POLYGON") %>%
  select("NAME_E") %>%
  mutate(area_km2 = as.numeric(st_area(geometry))/1000000) %>%
  filter(area_km2 > 1)

reserves_fda$NAME_E <- c("Tombstone_1", "Tombstone_2", "Tombstone_3")
names(reserves_fda)[names(reserves_fda) == "NAME_E"] <- "reserve"
reserves_fda <- reserves_fda[c("reserve", "area_km2", "geometry")]
reserves_sample <- reserves_fda
row.names(reserves_sample) <- c("1","2","3")

saveRDS(reserves_sample, file = "data/reserves_sample.rds")


# get streams in FDA 09EA
streams <- sf::st_read("C:/Users/mehou10/BEACONs Dropbox/Mélina Houle/gisdata/hydrology/borealC_v1_network.shp")

streams_sample <- streams %>%
  st_intersection(fda) %>%
  dplyr::select(UID_V6, BASIN, STRAHLER, ORDER1, ORDER2, ORDER3, UID)

saveRDS(streams_sample, file = "data/streams_sample.rds")
