## Started 13 March 2024 ##
## Satrted by Lizzie cribbing off Janneke's MastingAnalysis.R code ##

## JHRL says ...
# checkout README! and ...
# I think what we want is a figure showing that we have a long time series for many species OR many locations. 
# So perhaps AGO05 (PSME, TSHE, ABAM, THPL) would be good, or ABAM in any location it is at?

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/grephon/Masting/")
} else if (length(grep("ailene", getwd()))>0) 
{setwd("boomboom")
}

## 
## Start by cleaning the seeds! 
##

seedsthru2017 <- read.csv("./data/MORA_cleanseeds_2009-2017.csv", header=TRUE) #seeds
colnamesneeded <- names(seedsthru2017)[1:9] # removing totconeseeds (inelegantly)

seedsall <- seedsthru2017[0,]
for (year in c(2009:2022)){
  filehere <- read.csv(paste("./data/rawdata/sortedseeds/SortedSeeds_", year, ".csv", sep=""))
  filehere[["year"]] <- filehere[["Year"]]
  filehere[["stand"]] <- filehere[["Stand"]]
  filehere[["species"]] <- filehere[["Species"]]
  filehere[["trap"]] <- filehere[["Trapno"]]
  filehere[["filledseeds"]] <- filehere[["Whole_seeds"]]
  filehere[["emptyseeds"]] <- filehere[["Unfilled_seeds"]]
  filehere[["cones"]] <- filehere[["Cones"]]
  filehere[["conefilledseeds"]] <- filehere[["WholeSeeds_Cones"]]
  filehere[["coneemptyseeds"]] <- filehere[["UnfilledSeeds_Cones"]]
  if(year>2017){
    # Edits for 2018...
    filehere[["filledseeds"]] <- filehere[["Filled_seeds"]]
    filehere[["conefilledseeds"]] <- filehere[["FilledSeeds_Cones"]]
    }
  if(year>2020){
    # Edits for 2021...
    filehere[["species"]] <- filehere[["Spp"]]
    filehere[["filledseeds"]] <- filehere[["FilledSds"]] 
    filehere[["emptyseeds"]] <- filehere[["UnfilledSds"]]
    filehere[["conefilledseeds"]] <- filehere[["FilledSds.cones"]]
    filehere[["coneemptyseeds"]] <- filehere[["UnfilledSds.cones"]]
  }
  fileuse <- subset(filehere, select=colnamesneeded)
  seedsall <- rbind(seedsall, fileuse)
}

# Comparing suggests I need to add in ZEROS! 
if(FALSE){
subset(seedsthru2017, stand=="AG05" & year=="2009")
subset(seedsall, stand=="AG05" & year=="2009") 
}
# For now, do it cheaply using seedsthru2017 ... quite sure Janneke has a better way
seedsforzeroes <- expand.grid(year=c(2009:2022), stand=unique(seedsthru2017$stand), species=unique(seedsthru2017$species))
seedsallzeros  <- merge(seedsforzeroes, seedsall, all.x=TRUE, all.y=TRUE)
seedsallzeros[is.na(seedsallzeros)] <- 0
unique(seedsall$species) # Hmmm... I ignore this NA problem and error for now!
seedsallzeros$stand[seedsallzeros$stand=="PARA "] <- "PARA"


## 
## Now clean the traps 
## 

trapsize <- read.csv("./data/MORA_seedtrapinfo_2009-2017.csv", header=TRUE) #trapsizes

###reconfigure trapsize data
trapsize2 <- trapsize[,1:5]
size <- rep(NA, times=dim(trapsize)[1]); year <- size
trapsize2 <- cbind(trapsize[,1:5], year, size)
finaltrap <- c()

for(i in 2009:2017){
  trapsize3 <- trapsize2
  trapsize3[,6] <- rep(i, dim(trapsize2)[1])
  tmpsize <- trapsize[,i-2002]
  trapsize3[,7] <- tmpsize
  finaltrap <- rbind(finaltrap,trapsize3)
}

trapsizemerge <- trapsize

for (year in c(2018:2022)){
  filehere <- read.csv(paste("./data/rawdata/trapinfo/TrapSize_", year, ".csv", sep=""))
  filehere[["side"]] <- filehere[["Side"]]
  filehere[["year"]] <- filehere[["Year"]]
  filehere[["stand"]] <- filehere[["Stand"]]
  filehere[["trap"]] <- filehere[["Trap"]]
  filehere[["firstyeardata"]] <- filehere[["FirstYearData"]]
  fileheretomerge <- subset(filehere, select=c("side", "stand", "trap", "firstyeardata", "Size"))
  names(fileheretomerge)[names(fileheretomerge)=="Size"] <- paste("X", year, sep="")
  trapsizemerge <- merge(trapsizemerge, fileheretomerge, all.x=TRUE, all.Y=TRUE)
}

## ALERT! 
## Below is copied from Janneke's script and I don't actually know what it's doing...  
## reconfigure trapsize data
trapsize <- trapsizemerge
trapsize2 <- trapsize[,1:5]
size <- rep(NA, times=dim(trapsize)[1]); year <- size
trapsize2 <- cbind(trapsize[,1:5], year, size)
finaltrap <- c()

for(i in 2009:2022){
  trapsize3 <- trapsize2
  trapsize3[,6] <- rep(i, dim(trapsize2)[1])
  tmpsize <- trapsize[,i-2002]
  trapsize3[,7] <- tmpsize
  finaltrap <- rbind(finaltrap,trapsize3)
}

# Merge size and seed data
d <- merge(finaltrap, seedsallzeros, by = c("stand","trap","year"))
d$totfilledseeds <- d$filledseeds + d$conefilledseeds

## Plotting
##
library(ggplot2)

# Quick look by ABAM
dabam <- subset(d, species=="ABAM")
dabamagg <- aggregate(dabam[c("totfilledseeds")], dabam[c("stand", "year")], FUN=sum)

ggplot(dabamagg, aes(y=totfilledseeds, x=year, color=stand)) +
  geom_line() 

ggplot(dabamagg, aes(y=totfilledseeds, x=year)) +
  geom_line() + 
  facet_wrap(stand~.)

# Quick look for AGO5
specieshere <- c("ABAM", "PSME", "TSHE", "THPL")
dag05 <- subset(d[which(d$species %in% specieshere),], stand=="AG05")
dag05agg <- aggregate(dag05[c("totfilledseeds")], dag05[c("stand", "year", "species")], FUN=sum)

ggplot(dag05agg, aes(y=totfilledseeds, x=year, color=species)) +
  geom_line() 