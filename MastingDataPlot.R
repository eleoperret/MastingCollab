## Started 13 March 2024 ##
## Satrted by Lizzie who made some VERY quick plots ##

## Search for: Alert! in MastingDataMerge.R and make sure they are okay


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

# Quick look by ABAM
dabam <- subset(d, species=="ABAM")
dabamagg <- aggregate(dabam[c("totfilledseeds")], dabam[c("stand", "year")], FUN=sum)
dabamaggall <- aggregate(dabam[c("allseeds")], dabam[c("stand", "year")], FUN=sum)


ggplot(dabamagg, aes(y=totfilledseeds, x=year, color=stand)) +
  geom_line() + 
  scale_colour_viridis_d() + 
  theme_bw()

abamplotfacetfill <- ggplot(dabamagg, aes(y=totfilledseeds, x=year)) +
  geom_line() + 
  facet_wrap(stand~.) + 
  theme_bw()

abamplotfacetall <- ggplot(dabamaggall, aes(y=allseeds, x=year)) +
  geom_line() + 
  facet_wrap(stand~.) + 
  theme_bw()

# Quick look for AGO5
specieshere <- c("ABAM", "PSME", "TSHE", "THPL")
dag05 <- subset(d[which(d$species %in% specieshere),], stand=="AG05")
dag05agg <- aggregate(dag05[c("totfilledseeds")], dag05[c("stand", "year", "species")], FUN=sum)
dag05aggall <- aggregate(dag05[c("allseeds")], dag05[c("stand", "year", "species")], FUN=sum)


ag054sppfill <- ggplot(dag05agg, aes(y=totfilledseeds, x=year, color=species)) +
  geom_line() + 
  scale_colour_viridis_d(option="turbo") + 
  theme_bw()

ag054sppall <- ggplot(dag05aggall, aes(y=allseeds, x=year, color=species)) +
  geom_line() + 
  scale_colour_viridis_d(option="turbo") + 
  theme_bw()

ggsave("./figures/ag054sppfilled.pdf", ag054sppfill, height=4, width=6)
ggsave("./figures/ag054sppall.pdf", ag054sppall, height=4, width=6)
ggsave("./figures/abamplotfacetall.pdf", abamplotfacetall, height=6, width=8)
ggsave("./figures/abamplotfacetfill.pdf", abamplotfacetfill, height=6, width=8)