# This code is created to fit basic distributions
# Authors: Andreas K. Salk & Frederik Toftegaard

# Function to Install and Load R Packages
inst_load <- function(Required_Packages){
  Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
  
  if(length(Remaining_Packages)) 
  {
    install.packages(Remaining_Packages);
  }
  for(package_name in Required_Packages)
  {
    library(package_name,character.only=TRUE,quietly=TRUE);
  }
}

# The Packages needed to run the code 
Required_Packages=c("openxlsx", "fBasics", "MASS", "data.table", "triangle", "shinydashboard", "shiny", "rstudioapi");
# Call the Function
inst_load(Required_Packages)

# Changing working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# This will load the bestfit program ccreated to fit distributions and a function to check for package requirements
source("src/bestfit.R")

#Loading the data into R (You have to change the filename to fit your own placement of the file)
dat <- read.xlsx("Input-R.xlsx" ,sheet = 1, colNames = TRUE)

#Extracting the second column of the data
dat1 = dat[,2]

bestfit(dat1,15)