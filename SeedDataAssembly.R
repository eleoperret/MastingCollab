########################################################################
#### MORA seed data - pull together info needed for modeling        ####
#### Includes sorted seed data, seed trap size data, germinant data ####
#### Adds zeroes, and indicates where data is missing (NAs)         ####
########################################################################

###Files and paths to read in data
#seeds
seed_files <- list.files(path="./data/rawdata/sortedseeds/clean&notes")
rm_fls <- grep("error", seed_files) #remove error files
seed_files <- seed_files[-rm_fls]
seed_path <- "./data/rawdata/sortedseeds/clean&notes/"

#germinants in traps
germ_files <- list.files(path="./data/rawdata/germinantspertrap/clean&notes")
rm_flg <- grep("error", germ_files) #remove error files
germ_files <- germ_files[-rm_flg]
germ_path <- "./data/rawdata/germinantspertrap/clean&notes/"

#trap size, missing traps
info_files <- list.files(path="./data/rawdata/trapinfo")
ts_files <- grep("TrapSize", info_files) #identify trap sizes
trapsize_files <- info_files[ts_files]
inv_files <- grep("StandTrap", info_files) #identify inventory files
inventory_files <- info_files[inv_files]
info_path <- "./data/rawdata/trapinfo/"

###COMBINE DATA FILES
#Seeds
SortedSeeds_all <- c()

for(i in 1:length(seed_files)){
  tmpseed <- read.csv(paste(seed_path,seed_files[i],sep=""), header=TRUE)
  SortedSeeds_all <- rbind(SortedSeeds_all, tmpseed)
}

#Combine all germinant files
TrapGerms_all <- c()

for(i in 1:length(germ_files)){
  tmpgerm <- read.csv(paste(germ_path,germ_files[i],sep=""), header=TRUE)
  TrapGerms_all <- rbind(TrapGerms_all, tmpgerm)
}

#Combine Trap info
TrapSize_all <- c()

for(i in 1:length(trapsize_files)){
  tmpsize <- read.csv(paste(info_path,trapsize_files[i],sep=""), header=TRUE)
  yr0 <- strsplit(trapsize_files[i],"_")[[1]][2]
  yr <- as.numeric(substr(yr0,1,4))
  yrsize <- paste("X", yr, sep="")
  year <- rep(yr, length=dim(tmpsize)[1])
  tmpsize <- cbind(year, tmpsize)
  hdr <- dimnames(tmpsize)[[2]]
  hdr2 <- hdr
  hdr2[hdr=="Stand"]<-"stand"
  hdr2[hdr=="Trap"]<-"trapno"
  hdr2[hdr==yrsize]<-"size"
  hdr2[hdr=="Size"]<-"size"
  dimnames(tmpsize)[[2]] <- c(hdr2)
  #rm dside, firstyear, lastyeardata
  tmpsize <- tmpsize[,-(which(dimnames(tmpsize)[[2]]=="Side"))]
  tmpsize <- tmpsize[,-(which(dimnames(tmpsize)[[2]]=="FirstYearData"))]
  tmpsize <- tmpsize[,-(which(dimnames(tmpsize)[[2]]=="LastYearData"))]
  
  #now remove NAs
  tmpsize <- tmpsize[is.na(tmpsize$size)==FALSE,]
  
  #print(trapsize_files[i]); print(head(tmpsize))
  
  #now append
  TrapSize_all <- rbind(TrapSize_all, tmpsize)
}

#check Trapsize file DONE
#write.csv(TrapSize_all, file = paste(info_path,"Checkall.csv", sep=""), 
#          row.names=FALSE, col.names=TRUE)

#Remove traps that were lost, not collected
for(i in 1:length(inventory_files)){
  tmpinv <- read.csv2(paste(info_path,inventory_files[i],sep=""), sep=";", header=TRUE)
  allgood <- grep("no missing trap", tmpinv[1,1]) #check if no missing
  if(length(allgood)==1){next} #go to next if no missing traps
  
  for(j in 1:dim(tmpinv)[1]){
    tmpstnd <- tmpinv$stand[j]
    tmptrp <- tmpinv$trapno[j]
    tmpyr <- tmpinv$year[j]
    toremove <- which(TrapSize_all$stand==tmpstnd
                      &TrapSize_all$trapno==tmptrp
                      &TrapSize_all$year==tmpyr)
    TrapSize_all <- TrapSize_all[-toremove,]
  }
}
  

  
####Now assemble master file with seed data
#using data frames SortedSeeds_all, TrapGerms_all, TrapSize_all
#includes all years and traps sorted (minus missing traps)
#adds zeroes where trap sorted but no seeds found
#includes all species (regardless of whether they are in the stand or not)
#file is called SeedData_all

SeedData_all <- matrix(NA, nrow=0, ncol=12)
hdr1 <- c(dimnames(TrapSize_all)[[2]], dimnames(SortedSeeds_all)[[2]][4:10],
           dimnames(TrapGerms_all)[[2]][5])
dimnames(SeedData_all)[[2]] <- hdr1

##Now create a for loop that appends to SeedData_all
#Use TrapSize_all as a template for each species
#go through the years to fill non zero entries

allspp <- unique(SortedSeeds_all$spp) #first identify unique species, can change

for(i in 1:length(allspp)){
  tmpsize <- TrapSize_all
  spid <- rep(allspp[i], times=dim(tmpsize)[1])
  zeroes <- matrix(0, nrow=dim(tmpsize)[1],ncol=7)
  tmpsppdat <- cbind(tmpsize, spid, zeroes)
  dimnames(tmpsppdat)[[2]] <- hdr1
  
  #replace zeroes with NAs for cone seeds filled and unfilled from 2009-2012
  tmpsppdat$cones_sds_filled <- replace(tmpsppdat$cones_sds_filled, tmpsppdat$year<2013,NA)
  tmpsppdat$cones_sds_unfilled <- replace(tmpsppdat$cones_sds_filled, tmpsppdat$year<2013,NA)
  
  #extract species specific seed data, germination data
  tmpseed <- SortedSeeds_all[SortedSeeds_all$spp==allspp[i],]

  #for loop - goes through each row of tmpseed, ids tmpsppdat row to add info
  for(j in 1:dim(tmpseed)[1]){
    rwadd <- which(tmpsppdat$year==tmpseed$year[j] &
                   tmpsppdat$stand==tmpseed$stand[j] &
                   tmpsppdat$trapno==tmpseed$trapno[j])
    tmpsppdat[rwadd,6:11] <- tmpseed[j,5:10]  
  }
  
  #extract species specific seed data, germination data
  tmpgerm <- TrapGerms_all[TrapGerms_all$spp==allspp[i],]
  
  #for loop - goes through each row of tmpgerm, ids tmpsppdat row to add info
  for(j in 1:dim(tmpgerm)[1]){
    rwadd <- which(tmpsppdat$year==tmpgerm$year[j] &
                     tmpsppdat$stand==tmpgerm$stand[j] &
                     tmpsppdat$trapno==tmpgerm$trapno[j])
    tmpsppdat[rwadd,12] <- tmpgerm[j,5]  
  }
  
  #now append
  SeedData_all <- rbind(SeedData_all, tmpsppdat)
}

##Sorted seeds plus germinants
total_viable_sds <- SeedData_all$loose_sds_filled + SeedData_all$germintrap
mean(total_viable_sds)
SeedData_all <- cbind(SeedData_all, total_viable_sds)

#Write data


#Create year and stand and species specific means
SeedMeans <- tapply(total_viable_sds, list(as.factor(SeedData_all$spp), 
                                        as.factor(SeedData_all$stand),
                                        as.factor(SeedData_all$year)), FUN=mean)

#make graphs
par(mfrow=c(2,1), omi=c(0.2,0.2,0.2,0.2), mai=c(0.4,0.4,0.4,0.3), 
    mgp=c(1.25,0.3,0), xpd=FALSE, tck=-0.01)

#Abam
Abamseeds <- SeedMeans[1,,] + 0.08
mx <- max(Abamseeds, na.rm=TRUE)*1.05
yrs <- as.numeric(dimnames(Abamseeds)[[2]])

plot(yrs, Abamseeds[1,], xlim=c(2008.5,2024.5), ylim=c(0.08, mx), type="b",
     pch=21, bg="lightgreen", log="y", ylab="Log(seeds)", xlab="")

for(i in 2:dim(Abamseeds)[1]){
  lines(yrs, Abamseeds[i,])  
  points(yrs, Abamseeds[i,], pch=21, bg="lightgreen")
}
title("ABAM")

#Tshe
Tsheseeds <- SeedMeans[16,,] + 0.08
mx <- max(Tsheseeds, na.rm=TRUE)*1.05
yrs <- as.numeric(dimnames(Tsheseeds)[[2]])

plot(yrs, Tsheseeds[1,], xlim=c(2008.5,2024.5), ylim=c(0.08, mx), type="b",
     pch=21, bg="tan", log="y", ylab="Log(seeds)", xlab="")

for(i in 2:dim(Tsheseeds)[1]){
  lines(yrs, Tsheseeds[i,])  
  points(yrs, Tsheseeds[i,], pch=21, bg="tan")
}
title("TSHE")

##Stand Graphs - TO04, AM16
#
par(mfrow=c(2,1), omi=c(0.2,0.2,0.2,0.2), mai=c(0.4,0.4,0.4,0.3), 
    mgp=c(1.25,0.3,0), xpd=FALSE, tck=-0.01)

#TO04
TO04seeds <- SeedMeans[c(1,13,15,16),17,] + 0.08
mx <- max(TO04seeds, na.rm=TRUE)*1.5
yrs <- as.numeric(dimnames(TO04seeds)[[2]])

plot(yrs, TO04seeds[1,], xlim=c(2008.5,2024.5), ylim=c(0.08, mx), type="b",
     pch=21, bg="lightgreen", log="y", ylab="Log(seeds)", xlab="")

pltcols <- c("lightgreen", "darkgreen", "orange", "tan")

for(i in 2:dim(TO04seeds)[1]){
  lines(yrs, TO04seeds[i,])  
  points(yrs, TO04seeds[i,], pch=21, bg=pltcols[i])
}
title("TO04")
legend(x="topleft", c("ABAM", "PSME", "THPL", "TSHE"), 
       cex=0.65, pch=21, pt.bg=pltcols, horiz=TRUE, pt.cex = 1.15)

#AM16
AM16seeds <- SeedMeans[c(1,8,16,17),4,] + 0.08
mx <- max(AM16seeds, na.rm=TRUE)*2.5
yrs <- as.numeric(dimnames(AM16seeds)[[2]])

plot(yrs, AM16seeds[1,], xlim=c(2008.5,2024.5), ylim=c(0.08, mx), type="b",
     pch=21, bg="lightgreen", log="y", ylab="Log(seeds)", xlab="")

pltcols <- c("lightgreen", "lightblue", "tan", "plum")

for(i in 2:dim(AM16seeds)[1]){
  lines(yrs, AM16seeds[i,])  
  points(yrs, AM16seeds[i,], pch=21, bg=pltcols[i])
}
title("AM16")
legend(x="topleft", c("ABAM", "CANO", "TSHE", "TSME"), 
       cex=0.65, pch=21, pt.bg=pltcols, horiz=TRUE, pt.cex = 1.15) 


###To Do
#NAs for cone filled and unfilled (years where not separated)
#read out data

