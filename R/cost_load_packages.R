# RMRS 
# Authors: Jamie Sanderlin 
# Date: 1/26/2026

# Load packages (check if installed and if not install them)

################################################################################
## Load required packages
# Function from https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## First specify the packages of interest
packages = c("here","dplyr", "lubridate","tidyverse", 'readxl',
             "reshape2","ggplot2","psych","wesanderson")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)