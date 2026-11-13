# reserve is added as expected

    Code
      append_reserve(benchmarks, "network", pa_1, "PA_1")
    Output
      Simple feature collection with 3 features and 1 field
      Geometry type: POLYGON
      Dimension:     XY
      Bounding box:  xmin: 478494.8 ymin: -1540493 xmax: 752950.9 ymax: -1305868
      Projected CRS: NAD_1983_Albers
        network                       geometry
      1 PB_0001 POLYGON ((689337.4 -1445340...
      2 PB_0002 POLYGON ((686200.9 -1499500...
      3    PA_1 POLYGON ((478494.8 -1382608...

# reserve is dissolved correctly

    Code
      append_reserve(benchmarks, "network", pa_1, "PA_1")
    Output
      Simple feature collection with 2 features and 1 field
      Geometry type: POLYGON
      Dimension:     XY
      Bounding box:  xmin: 686200.9 ymin: -1540493 xmax: 752950.9 ymax: -1418750
      Projected CRS: NAD_1983_Albers
        network                       geometry
      1 PB_0001 POLYGON ((689337.4 -1445340...
      2    PA_1 POLYGON ((686201 -1499500, ...

# networks made as expected

    Code
      build_network_polygons(benchmarks, "network", c("PB_0001__PB_0002",
        "PB_0001__PB_0003"))
    Output
      Simple feature collection with 2 features and 1 field
      Geometry type: GEOMETRY
      Dimension:     XY
      Bounding box:  xmin: 686200.9 ymin: -1540493 xmax: 752950.9 ymax: -1418750
      Projected CRS: NAD_1983_Albers
                 network                       geometry
      1 PB_0001__PB_0002 MULTIPOLYGON (((686201 -149...
      2 PB_0001__PB_0003 POLYGON ((715701.1 -1487750...

# upstream catchments dissolve correctly

          PA    AWI                       geometry
    1 191337 0.5000 POLYGON ((703200.9 -1476500...
    2 191396 0.7399 POLYGON ((726950.9 -1473750...
    3 199392 0.8446 POLYGON ((681700.9 -1490750...

# benchmarks build as expected

    Code
      dissolve_catchments_from_table(catchments_sample, builder_table_sample,
        "network", dissolve_list = c("PB_0001", "PB_0002"))
    Output
      Simple feature collection with 2 features and 1 field
      Geometry type: POLYGON
      Dimension:     XY
      Bounding box:  xmin: 686200.9 ymin: -1540493 xmax: 752950.9 ymax: -1418750
      Projected CRS: NAD_1983_Albers
        network                       geometry
      1 PB_0001 POLYGON ((689337.4 -1445340...
      2 PB_0002 POLYGON ((686200.9 -1499500...

# networks build as expected

    Code
      dissolve_catchments_from_table(catchments_sample, builder_table_sample,
        "network", dissolve_list = c("PB_0001__PB_0002"))
    Output
      Simple feature collection with 1 feature and 1 field
      Geometry type: MULTIPOLYGON
      Dimension:     XY
      Bounding box:  xmin: 686200.9 ymin: -1540493 xmax: 752950.9 ymax: -1418750
      Projected CRS: NAD_1983_Albers
                 network                       geometry
      1 PB_0001__PB_0002 MULTIPOLYGON (((689337.4 -1...

