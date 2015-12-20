# Run all examples

# Create the output directory
# Notice that this directory is defined again in all main_*.R files as OUTPUT.DIR
dir.create("Output", showWarnings = FALSE, recursive = TRUE)

source("main_inventory.R")
source("main_stocks.R")
source("main_stocks_compression_tests.R")



