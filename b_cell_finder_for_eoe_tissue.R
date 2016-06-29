#######################################################################
# R code to Identify B cells in EoE Tissue Plates 
# Created by Cecilia Washburn
# June 2016

# The EoE T cell profiling project single-cell-sorts plates of digested tissue
# The gate includes T cells and B cells
# The machine collects index data that can be used to identify cell phenotype after the fact

# This script reads the index sort .csv file for individual plates and prints wells that are CD19+/CD3- (B cells) and wells that are CD19+/CD3+ (possibly indicating a dead cell)
# The function also optionally prints a few graphs of the data

#Part 1 (line 26) defines function b_cell_finder and runs on one file with graphs
#Part 2 (line 124) prints the Bcell locations of a full directory of files to the consol
#Part 3 (line 135) contains a code to collate all .csv index sort files to find trends

#To Use:
##  Set working directory to appropriate file path (line 38)
##  Set 'csv' as file of interest (line 41)
##  Run lines 25-119

#### Script will need to be adjusted if a new antibody panel is used.

#######################################################################

################
#Part one: b_cell_finder function and individual file reading
###############

#Get some things set up
library(dplyr)
library(ggplot2)

rm(list=ls())

###Set working directory to where your files are stored
#OR include the full file path below
setwd("C:/Users/cew27/Documents/Index Sort")

#your file name goes here
csv <- "06-09-2016 plate 2 72.csv"

#Function
#Graphs argument determines if you would want graphs of the data printed. Default is True, but can be turned off to run a list of files faster.

b_cell_finder <- function(csv, graphs = T){
  
  #Reads file, pares down to CD19 and CD3 data, then converts to numeric
  plate1 <- read.csv(csv, skip = 15, header = T, dec = ".")

  plate1 <- select(plate1, Well, CD19 = All.Events.CD19.Alexa.Fluor.700.A.Mean,
                   CD3 = All.Events.CD3.FITC.A.Mean )

  plate1$CD3 <- gsub(",", "", plate1$CD3)
  plate1$CD19 <- gsub(",", "", plate1$CD19)
  
  plate1[2:3] <- sapply(plate1[2:3], as.numeric)


  #let's assume a negative number is negative for that marker and leave it at zero
  plate1$CD3[plate1$CD3 <0] <- 0
  plate1$CD19[plate1$CD19 <0] <- 0

  #Graph of CD19 activity
  if(graphs == T) {
    print(
      ggplot(plate1, aes(x = Well, y = CD19))+
        geom_bar(stat = "identity", fill = "skyblue1")+
        ggtitle("CD19 Expression by well")+
        ylab("CD19")
   )}
  
  #Graph of CD19 vs CD3
  if(graphs == T) {
    print(
      ggplot(plate1, aes(y = CD3, x = CD19)) +
        geom_jitter(color = "darkmagenta", fill= "darkmagenta", shape = 23)+
        scale_x_log10()+
        scale_y_log10()+
        ggtitle("Index Sorted Cells")+
        ylab("CD3")+
        xlab("CD19")
  ) }
  
  #Find the cells that were CD19+
  bcells <- filter(plate1, CD19 > 100 & CD3 < 1000)

  if(graphs == T) {
    print(
      ggplot(bcells, aes(y = CD3, x = CD19)) +
        geom_jitter(color = "darkmagenta", fill= "darkmagenta", shape = 23)+
        scale_x_log10()+
        scale_y_log10()+
        ggtitle("Bcells")+
        ylab("CD3")+
        xlab("CD19")
  )}
  
  #Find dead cells (double positives)
  dead <- filter(plate1, CD19 > 100 & CD3 >1000)
  
  if(graphs == T) {
    print(
      ggplot(dead, aes(y = CD3, x = CD19)) +
        geom_jitter(color = "black", fill= "black", shape = 23)+
        scale_x_log10()+
        scale_y_log10()+
        ggtitle("Dead cells")+
        ylab("CD3")+
        xlab("CD19")
  )}
  
  #print well numbers to the console
  AllB <- c("B cells: ", as.character(bcells$Well))
  print(paste(AllB, collapse = " "))
  AllDead <- c("'Dead' (double positive) cells: ", as.character(dead$Well))
  return(paste(AllDead, collapse = " "))
  
}

#Run the function on your file
b_cell_finder(csv, graphs = T)

#######################
#Part 2
#Code to run to print a list of b cell wells by file
#Assumes that your working directory is only index sort files
files <- list.files() 

for(i in 1:length(files)){
  print(files[i])
  print(b_cell_finder(files[i], graphs = F))}


########################################
#Part 3
### Visualizing index data as a whole
#This is where the cut-off for CD19+ and CD3+ comes from

#Lists the index files from a folder containing all tissue index files
files <- list.files() 
filenum <- 0

#merges all files into single data frame
for (i in 1:length(files)) {
  
  plate <- read.csv(files[i], skip = 15, header = T, dec = ".")
  plate <- select(plate, CD19 = All.Events.CD19.Alexa.Fluor.700.A.Mean, 
                  CD3 = All.Events.CD3.FITC.A.Mean)
  
  #Converts data to numeric
  plate$CD3 <- gsub(",", "", plate$CD3)
  plate$CD19 <- gsub(",", "", plate$CD19)
  
  plate[1:2] <- sapply(plate[1:2], as.numeric)

  if(filenum == 0){
    all <- plate
    filenum <- filenum + 1}
  
  else {all <- rbind(all, plate)}
}

#let's assume a negative number is negative for that marker and leave it at zero
all$CD3[all$CD3 <0] <- 0
all$CD19[all$CD19 <0] <- 0

#Graph of CD3 and CD19
ggplot(all, aes(y = CD3, x = CD19)) +
  geom_point(color = "darkmagenta", fill= "darkmagenta", shape = 23)+
  scale_x_log10()+
  scale_y_log10()+
  stat_density2d(contour = TRUE, n = 100)+
  ggtitle("All Index Sorted Cells")+
  ylab("CD3+")+
  xlab("CD19+")

#Data frames of all CD19+ and CD19+/CD3+ (dead) cells
totalb <- filter(all, CD19 > 100 & CD3 < 1000)

totaldead <- filter(all, CD19 > 100 & CD3 >1000)

