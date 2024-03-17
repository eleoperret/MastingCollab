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
library(viridis)

d <- read.csv("./output/seedswithtraps_allyears_quickdirty.csv")
standz <- read.csv("./data/MORA STANDS_LatLongs_UTM.csv")

# Quick look by ABAM
dabam <- subset(d, species=="ABAM")
dabamagg <- aggregate(dabam[c("totfilledseeds")], dabam[c("stand", "year")], FUN=sum)
dabamaggall <- aggregate(dabam[c("allseeds")], dabam[c("stand", "year")], FUN=sum)

## ggplot versions

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
specieshere <- c("ABAM", "PSME",  "THPL", "TSHE")
latbihere <- c(expression(italic("Abies amabilis")), 
    expression(italic("Pseudotsuga menziesii")), 
    expression(italic("Thuja plicata")),
    expression(italic("Tsuga heterophylla")))

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

## Base R versions 


## One stand
colz <- viridis(6)

pdf("./figures/ag054sppfilledbaseR.pdf", height=6, width=8)
par(mar=c(4,5,2,2))
plot(totfilledseeds~year, data=dag05agg, type="n", xlab="Year", 
  ylab=expression(paste("Total filled seeds (m"^"-2", ")", sep="")), ylim=c(0, 10800)) 
for (i in 1:length(specieshere)){
  subby <- subset(dag05agg, species==specieshere[i])
  lines(totfilledseeds~year, data=subby, col=colz[i], lwd=2)
  points(totfilledseeds~year, data=subby, col=colz[i], pch=16)
}
legend("topleft", latbihere, pch=16, col=colz[1:4], bty="n", lwd=2)
dev.off()

pdf("./figures/ag054sppfilledlnbaseR.pdf", height=6, width=8)
par(mar=c(4,5,2,2))
plot(log(totfilledseeds)~year, data=dag05agg, type="n", xlab="Year", 
  ylab=expression(paste("ln(Total filled seeds (m"^"-2", ")", sep="")), ylim=c(3,11)) 
for (i in 1:length(specieshere)){
  subby <- subset(dag05agg, species==specieshere[i])
  lines(log(totfilledseeds)~year, data=subby, col=colz[i], lwd=2)
  points(log(totfilledseeds)~year, data=subby, col=colz[i], pch=16)
}
legend("topleft", latbihere, pch=16, col=colz[1:4], bty="n", lwd=2)
dev.off()

## One species
colz <- rev(viridis(nrow(standz)))

pdf("./figures/abamplotbaseR.pdf", height=6, width=8)
par(mar=c(4,5,2,2))
plot(totfilledseeds~year, data=dabamagg, type="n", xlab="Year", 
  ylab=expression(paste("Total filled seeds (m"^"-2", ")", sep=""))) 
for (i in 1:nrow(standz)){
  subby <- subset(dabamagg, stand==standz$PlotNo[i])
  lines(totfilledseeds~year, data=subby, col=colz[i], lwd=2)
  points(totfilledseeds~year, data=subby, col=colz[i], pch=16)
}
legend("topleft", standz$PlotNo, pch=16, col=colz, bty="n", lwd=2)
dev.off()

pdf("./figures/abamplotlnbaseR.pdf", height=6, width=8)
par(mar=c(4,5,2,7), xpd=TRUE)
plot(log(totfilledseeds)~year, data=dabamagg, type="n", xlab="Year", 
  ylab=expression(paste("ln(Total filled seeds (m"^"-2", "))", sep=""))) 
for (i in 1:nrow(standz)){
  subby <- subset(dabamagg, stand==standz$PlotNo[i])
  lines(log(totfilledseeds)~year, data=subby, col=colz[i], lwd=2)
  points(log(totfilledseeds)~year, data=subby, col=colz[i], pch=16)
}
legend("topright", inset=c(-0.2,0), standz$PlotNo, pch=16, col=colz, bty="n", lwd=2)
dev.off()



standzsm <- standz[which(standz$PlotNo %in% c("AG05", "AO03", "AV06", "AV14", "PARA")),]
colz <- rev(viridis(nrow(standzsm)))

pdf("./figures/abamplotlnbaseR_fewerstands.pdf", height=6, width=8)
par(mar=c(4,5,2,7), xpd=TRUE)
plot(log(totfilledseeds)~year, data=dabamagg, type="n", xlab="Year", 
  ylab=expression(paste("ln(Total filled seeds (m"^"-2", "))", sep=""))) 
for (i in 1:nrow(standz)){
  subby <- subset(dabamagg, stand==standzsm$PlotNo[i])
  lines(log(totfilledseeds)~year, data=subby, col=colz[i], lwd=2)
  points(log(totfilledseeds)~year, data=subby, col=colz[i], pch=16)
}
legend("topright", standzsm$PlotNo, pch=16, col=colz, bty="n", lwd=2)
dev.off()