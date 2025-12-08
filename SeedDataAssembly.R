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

for(i in 1:length(seed_files)){
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



