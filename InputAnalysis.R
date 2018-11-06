# This code is created to fit basic distributions
# Authors: Andreas K. Salk & Frederik Toftegaard

# Set you working directory (File path where the InputAnalyzer folder is saved)
setwd("C:/Users/andre/Downloads/InputAnalyzer/src")
# This will load the bestfit program ccreated to fit distributions and a function to check for package requirements
source("bestfit.R"); source("Install_And_Load.R")
# The Packages needed to run the code 
Required_Packages=c("openxlsx", "fBasics", "MASS", "data.table", "triangle", "shinydashboard", "shiny");
# Call the Function
Install_And_Load(Required_Packages)

#Loading the data into R (You have to change the filename to fit your own placement of the file)
dat <- read.xlsx("Input-R.xlsx" ,sheet = 1, colNames = TRUE)

#Extracting the second column of the data
dat1 = dat[,2]

bestfit(dat1,15)

