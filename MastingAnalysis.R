##############################################################
#### Script to analyze seed masting data from Mt. Rainier ####
#### Script created 12/01/2024; last worked on 12/01/2024 ####
##############################################################

###Set directory
setwd("C:/Users/janne/Dropbox/Research/WA_Projects/MtRainier/seeds&seedlings/2024_Data Assembly/Analysis")

###First read in seed data and seed trap size data
###Note - seed data have been cleaned and compiled from individual years
seeds <- read.csv("MORA_cleanseeds_2009-2017.csv", header=TRUE) #seeds
trapsize <- read.csv("MORA_seedtrapinfo_2009-2017.csv", header=TRUE) #trapsizes

