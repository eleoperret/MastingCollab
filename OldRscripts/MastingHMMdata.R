## Started 12 Mar 2025 ##
## By Lizzie ##

## Grab data for Mike ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/grephon/Masting/")
} else if (length(grep("ailene", getwd()))>0) 
{setwd("boomboom")
}

# libraries
library(ggplot2)

d <- read.csv("./output/seedswithtraps_allyears_quickdirty.csv")

# Subset to TSHE
tshe <- subset(d, species=="TSHE")

ggplot(tshe, aes(filledseeds)) +
	geom_histogram() +
	facet_wrap(.~stand)

# Based on my notes in MastingDataMerge.R and the above plot, I think these data have the zeros...
# But I am not sure about TSHE
mospplist <- c("TSHE", "THPL", "ABAM", "PSME")
mospp <- d[which(d$species %in% mospplist),]

ggplot(mospp, aes(filledseeds)) +
	geom_histogram() +
	facet_grid(species~stand)

# Okay, will stick with TSHE
tsheav06 <- subset(tshe, stand=="AV06")
tsheav06$allseeds <- tsheav06$filledseeds + tsheav06$emptyseeds
tsheav06sm <- subset(tsheav06, select=c("trap", "year", "size","allseeds"))


write.csv(tsheav06sm, "output/quickdatTSHEforMike.csv", row.names=FALSE)