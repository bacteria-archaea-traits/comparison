#################
# Load packages #
#################

library("ggpubr")
library("png")
library("readxl")
library("tidyverse")
library("Hmisc")
library("agricolae") #For New Tukeys method
library("caTools") #For regression analysis

options(stringsAsFactors = FALSE)

# Load custom functions for PCA
source("R/PCA_functions.R")

###########
# Process #
###########

# Prepare data frames
source("R/prep.R")

# Create figures
source("R/figures.R")
