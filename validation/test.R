# Script to validate disturbance attributes

# Load function
source('validate_functions.R')

# Function arguments

## data_pkg = name of data package
## lyr_line = name of linear disturbances layer
## lyr_poly = name of areal disturbances layer
## output_dir = location of output directory (will be created if it doesn't exist)

# Examples using FDA_10AA_001
validate_project(data_pkg='Data_package.gpkg', lyr_line='YG_Line', lyr_poly='YG_Poly', output_dir='output/')
