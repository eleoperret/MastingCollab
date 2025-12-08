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
Seeds_all <- c()

for(i in 1:length(seed_files)){
  tmpseed <- read.csv(paste(seed_path,seed_files[i],sep=""), header=TRUE)
  Seeds_all <- rbind(Seeds_all, tmpseed)
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
  
  print(trapsize_files[i]); print(head(tmpsize))
  
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
  }
  

  
  
  
  ####
##Determine unique species
allspp <- unique(Seeds_all$spp)