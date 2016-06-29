#######################################################################
# R code to Identify CRTH2+ T cells in EoE Overnight Culture Plates
# Created by Cecilia Washburn
# June 2016

# The EoE T cell profiling project single-cell-sorts plates of T cells that have been cultured in allergen overnight
# Cells are sorted on CD154 and CD69 expression
# The machine collects index data that can be used to identify CRTH2+ cells after the fact

# This script reads the index sort .csv file for individual plates and prints wells that express CRTH2
#The CRTH2 gate is variable: users have the option to set the cut-off for CRTH2 activity based on flowjo graphs of the CRTH2 FMO
# The function also optionally prints a few graphs of the data

#Part 1 (line 29) defines function crth2_finder and runs on one file with graphs
#Part 2 (line 124) prints the CRTH2+ wells of a full directory of files to the consol
#Part 3 (line 135) contains a code to collate all .csv index sort files to find trends

#To Use:
##  Set working directory to appropriate file path (line 39)
##  Set 'csv' as file of interest (line 42)
##  Set your gate limit from FMO data on FlowJo (line 45)
##  Run lines 32-122

#### Script will need to be adjusted if a new antibody panel is used.
#######################################################################


################
#Part one: crth2_finder function and individual file reading
###############

library(dplyr)
library(ggplot2)

rm(list=ls())

###Set working directory to where your files are stored
#OR include the full file path below
setwd("C:/Users/cew27/Documents/CRTH2 Index")

#your file name goes here
csv <- "EOE culture 02-12-2016 plate 247 2.csv"

#your CRTH2 gate limit goes here (use the CRTH2 FMO data from flowjo)
cr_fmo <- 100

#Function
#Graphs argument determines if you would want graphs of the data printed. Default is True, but can be turned off to run a list of files faster.

crth2_finder <- function(csv, cr_limit = 100, graphs = T){
  
  #Reads file, pares down to CD69, CD154, CD3 and CRTH2 data, then converts to numeric
  plate1 <- read.csv(csv, skip = 15, header = T, dec = ".")
  
  plate1 <- select(plate1, Well, CD69 = All.Events.CD69.BV711.A.Mean, 
                   CD154 = All.Events.CD154.APC.A.Mean, 
                   CRTH2 = All.Events.CRTH2.PE.Texas.Red.A.Mean,
                   CD3 = All.Events.CD3.BV650.A.Mean)
  
  plate1$CD69 <- gsub(",", "", plate1$CD69)
  plate1$CD154 <- gsub(",", "", plate1$CD154)
  plate1$CRTH2 <- gsub(",", "", plate1$CRTH2)
  plate1$CD3 <- gsub(",", "", plate1$CD3)
  
  plate1[2:5] <- sapply(plate1[2:5], as.numeric)
  
  
  #let's assume a negative number is negative for that marker and leave it at zero
  plate1$CD69[plate1$CD69 <0] <- 0
  plate1$CD154[plate1$CD154 <0] <- 0
  plate1$CRTH2[plate1$CRTH2 <0] <- 0
  plate1$CD3[plate1$CD3 <0] <- 0
  
  #Bar Graph of CRTH2 activity by well
  if(graphs == T) {
    print(
      ggplot(plate1, aes(x = Well, y = CRTH2))+
        geom_bar(stat = "identity", fill = "darkred")+
        ggtitle("CRTH2 Expression by well")
    )}
  
  #Graph of CD69 vs CD154
  if(graphs == T) {
    print(
      ggplot(plate1, aes(x = CD69, y = CD154)) +
        geom_jitter(color = "blue3", fill= "blue3", shape = 19)+
        scale_x_log10(limits = c(1, 100000))+
        scale_y_log10(limits = c(1, 100000))+
        stat_density2d(contour = TRUE, n = 100, na.rm = T)+
        ggtitle("CD69 vs CD154")
    ) }
  
  #Graph of CRTH2 activity by CD3
  if(graphs == T) {
    print(
      ggplot(plate1, aes(x = CRTH2, y = CD3)) +
        geom_jitter(color = "darkred", fill= "darkred", shape = 15)+
        scale_x_log10()+
        scale_y_log10()+
        ggtitle("CRTH2 Activity")
    ) }
  
  #Find the cells had CRTH2 activity
  crth2pos <- filter(plate1, CRTH2 > cr_limit)
  
   #print positive well numbers to the console
  active <- c("CRTH2 Positive: ", as.character(crth2pos$Well))
  print(length(active)-1)
  return(paste(active, collapse = " "))
  
}

#Run the function on your file
crth2_finder(csv, graphs = T, cr_limit = cr_fmo)

#######################
#Part 2
#Code to run to print a list of CRTH2+ wells by file
#Assumes that your working directory is only index sort files
files <- list.files() 

for(i in 1:length(files)){
  print(files[i])
  print(crth2_finder(files[i], graphs = F))}


########################################
#Part 3
### Visualizing index data as a whole

#Lists the index files from a folder containing all tissue index files
files <- list.files() 
filenum <- 0

#merges all files into single data frame
for (i in 1:length(files)) {
  
  plate <- read.csv(files[i], skip = 15, header = T, dec = ".")
  plate <- select(plate, CD69 = All.Events.CD69.BV711.A.Mean, 
                   CD154 = All.Events.CD154.APC.A.Mean, 
                   CRTH2 = All.Events.CRTH2.PE.Texas.Red.A.Mean,
                   CD3 = All.Events.CD3.BV650.A.Mean)
  
  plate$CD69 <- gsub(",", "", plate$CD69)
  plate$CD154 <- gsub(",", "", plate$CD154)
  plate$CRTH2 <- gsub(",", "", plate$CRTH2)
  plate$CD3 <- gsub(",", "", plate$CD3)
  
  plate[1:4] <- sapply(plate[1:4], as.numeric)
  
  if(filenum == 0){
    all <- plate
    filenum <- filenum + 1}
  
  else {all <- rbind(all, plate)}
}

#let's assume a negative number is negative for that marker and leave it at zero
all$CD69[all$CD69 <0] <- 0
all$CD154[all$CD154 <0] <- 0
all$CRTH2[all$CRTH2 <0] <- 0
all$CD3[all$CD3 <0] <- 0

#Graph of CD69 vs CD154
    ggplot(all, aes(x = CD69, y = CD154)) +
      geom_jitter(color = "blue3", fill= "blue3", shape = 19)+
      scale_x_log10()+
      scale_y_log10()+
      stat_density2d(contour = TRUE, n = 100) +
      ggtitle("CD69 vs CD154")

#Graph of CRTH2 activity by CD3
    ggplot(all, aes(x = CRTH2, y = CD3)) +
      geom_jitter(color = "darkred", fill= "darkred", shape = 15)+
      scale_x_log10()+
      scale_y_log10()+
      ggtitle("CRTH2 Activity")

#Data frame of all CRTH2+ cells
totalb <- filter(all, CD19 > 100 & CD3 < 1000)


