# testthat runs this file before the tests.

# Many tests call sf functions (st_transform, st_area, st_crs, ...) without sf::
library(sf)

# Sample objects used by the tests: catchments_sample, builder_table_sample,
# led_sample, streams_sample, habitat_sample, ref_poly
# Helpers run from tests/testthat, both in devtools::test() and devtools::load_all()
for (f in list.files("testdata", "[.]rda$", full.names = TRUE)) load(f)
led_sample     <- terra::rast("testdata/led_sample.tif")
habitat_sample <- terra::rast("testdata/habitat_sample.tif")

rm(f)
