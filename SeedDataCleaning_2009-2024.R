####################################################################
### script to clean seed files from 2009-2024                    ###
### also used to generate 2009-2019 StandTrapInventory Files     ###
### created by Janneke Hille Ris Lambers, 17.10.2025             ###
####################################################################

###Read in 2009 - 2024 seed trap data
###Check species, missing traps, etc
###creates errorcheck file
###seperate text file provides info on solutions
## names of data files

##Sorted seed data
seed_files <- list.files(path="./data/rawdata/sortedseeds")
seed_files <- seed_files[-1]
seed_path <- "./data/rawdata/sortedseeds/"
outputfile <- paste(seed_path,"clean&notes/errorcheck_10.2025.txt", sep="")
write("October 17, 2025 ERROR CHECK", outputfile)

for(i in 1:length(seed_files)){
  tempseed <- read.csv(paste(seed_path, seed_files[i], sep=""))
  
  #extract year
  tmp1 <- strsplit(seed_files[i],".csv")[[1]]
  yr <- strsplit(tmp1,"MORA_")[[1]][2]
  yr <- as.numeric(yr)
  
  #Signal new file
  write("NEW FILE CHECK", file=outputfile, append=TRUE)
  write(seed_files[i], file=outputfile, append=TRUE)
  
  #check years
  write("years", file=outputfile, append=TRUE)
  if(yr<2020){write(unique(tempseed$Year), file=outputfile, append=TRUE)}
  if(yr>2019){write(unique(tempseed$year), file=outputfile, append=TRUE)}

  #check stand
  write("stands", file=outputfile, append=TRUE)
  if(yr<2020){write(unique(tempseed$Stand), file=outputfile, append=TRUE)}
  if(yr>2019){write(unique(tempseed$stand), file=outputfile, append=TRUE)}
  
  #check species
  write("species", file=outputfile, append=TRUE)
  if(yr<2020){write(unique(tempseed$Species), file=outputfile, append=TRUE)}
  if(yr>2019){write(unique(tempseed$spp), file=outputfile, append=TRUE)}
  
  ###Now check stand / traps present in that year
  if(yr<2020){Stands <- unique(tempseed$Stand)}
  if(yr>2019){Stands <- unique(tempseed$stand)}
  
  write("Check traps per stand", file=outputfile, append=TRUE)
  
  for(j in 1:length(Stands)){
    if(yr<2020){tempseed_stand <- tempseed[tempseed$Stand==Stands[j],]}
    if(yr>2019){tempseed_stand <- tempseed[tempseed$stand==Stands[j],]}
    
    write(Stands[j], file=outputfile, append=TRUE)
    write("traps present", file=outputfile, append=TRUE)
    if(yr<2020){write(unique(tempseed_stand$Trapno), file=outputfile, append=TRUE)}
    if(yr>2019){write(unique(tempseed_stand$trapno), file=outputfile, append=TRUE)}
  }
}

##Now make corrections to individual sorted seed data, write to clean folder
##

#############################
##2009 DATA
Raw2009data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2009.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2009data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2009data)[[2]] <- c(hdr2)
Raw2009data <- Raw2009data[,-(which(hdr2=="Date"))]

##now write out
write.csv(Raw2009data, file = paste(seed_path,"clean&notes/",
                             "SortedSeeds_MORA_2009.csv", sep=""), 
                              row.names=FALSE, col.names=TRUE)

#############################
##2010 DATA
Raw2010data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2010.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2010data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2010data)[[2]] <- c(hdr2)
Raw2010data <- Raw2010data[,-(which(hdr2=="Date"))]

#Remove row with MISSING in species
rmrow <- which(Raw2010data$spp=="MISSING") #remove row with MISSING as species
Raw2010data <- Raw2010data[-rmrow,]

##now write out
write.csv(Raw2010data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2010.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)

#############################
##2011 DATA
Raw2011data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2011.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2011data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2011data)[[2]] <- c(hdr2)
Raw2011data <- Raw2011data[,-(which(hdr2=="Date"))]

#Remove row with MISSING in species
rmrow <- which(Raw2011data$spp=="NOSEEDS") #remove row indicating 0 seeds
Raw2011data <- Raw2011data[-rmrow,]

##now write out
write.csv(Raw2011data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2011.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)

#############################
#2012 DATA
Raw2012data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2012.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2012data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2012data)[[2]] <- c(hdr2)
Raw2012data <- Raw2012data[,-(which(hdr2=="Date"))]

##now write out
write.csv(Raw2012data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2012.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)

#############################
##2013 DATA
Raw2013data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2013.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2013data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2013data)[[2]] <- c(hdr2)
Raw2013data <- Raw2013data[,-(which(hdr2=="Date"))]

#Remove row with MISSING in species
rmrow <- which(Raw2013data$spp=="MYST") #remove row with mystery seeds
Raw2013data <- Raw2013data[-rmrow,]
rmrow <- which(Raw2013data$spp=="NOSEEDS") #remove row indicating trap has 0 seeds
Raw2013data <- Raw2013data[-rmrow,]

##now write out
write.csv(Raw2013data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2013.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)

#############################
##2014 DATA
Raw2014data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2014.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2014data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2014data)[[2]] <- c(hdr2)
Raw2014data <- Raw2014data[,-(which(hdr2=="Date"))]

#Remove row with MISSING in species
rmrow <- which(Raw2014data$spp=="MYST") #remove row with mystery seeds
Raw2014data <- Raw2014data[-rmrow,]

##now write out
write.csv(Raw2014data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2014.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)

#############################
##2015 DATA
Raw2015data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2015.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2014data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2015data)[[2]] <- c(hdr2)
Raw2015data <- Raw2015data[,-(which(hdr2=="Date"))]

##now write out
write.csv(Raw2015data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2015.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)

#############################
##2016 DATA
Raw2016data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2016.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2016data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2016data)[[2]] <- c(hdr2)
Raw2016data <- Raw2016data[,-(which(hdr2=="Date"))]

#Remove row with MISSING in species
rmrow <- which(Raw2016data$spp=="MISSING") #remove where seeds not collected
Raw2016data <- Raw2016data[-rmrow,]
rmrow <- which(is.na(Raw2016data$trapno)==TRUE) #remove where trap no unknow
Raw2016data <- Raw2016data[-rmrow,]

##now write out
write.csv(Raw2016data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2016.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)


#############################
##2017 DATA
Raw2017data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2017.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2017data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2017data)[[2]] <- c(hdr2)
#rm date header
Raw2017data <- Raw2017data[,-(which(hdr2=="Date"))]

#Remove row with MISSING in species
rmrow <- which(Raw2017data$spp=="MISSING") #remove where seeds not collected
Raw2017data <- Raw2017data[-rmrow,]
rmrow <- which(Raw2017data$spp=="NOSEEDS") #remove indicator of zero seeds
Raw2017data <- Raw2017data[-rmrow,]

#Change ACGL Acer glabrum, to ACCI
Raw2017data$spp[Raw2017data$spp[]=="ACGL"] <- "ACCI"

##now write out
write.csv(Raw2017data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2017.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)

#############################
##2018 DATA
Raw2018data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2018.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2018data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2018data)[[2]] <- c(hdr2)

#Remove row with MISSING in species
rmrow <- which(is.na(Raw2018data$spp)==TRUE) #remove NA species- means 0 seeds
Raw2018data <- Raw2018data[-rmrow,]

##now write out
write.csv(Raw2018data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2018.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)

#############################
##2019 DATA
Raw2019data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2019.csv",sep=""), header=TRUE)
hdr <- dimnames(Raw2019data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Whole_seeds"]<-"loose_sds_filled"
hdr2[hdr=="Unfilled_seeds"]<-"loose_sds_unfilled"
hdr2[hdr=="Cones"]<-"cones"
hdr2[hdr=="WholeSeeds_Cones"]<-"cones_sds_filled"
hdr2[hdr=="UnfilledSeeds_Cones"]<-"cones_sds_unfilled"
hdr2[hdr=="SeedsInCones"]<-"cones_sds_all"
hdr2[hdr=="Notes"]<-"standardized_notes"
dimnames(Raw2019data)[[2]] <- c(hdr2)

#Change ALSI, Alnus sitchensis, to ALVI
Raw2019data$spp[Raw2019data$spp[]=="ALSI"] <- "ALVI"

##now write out
write.csv(Raw2019data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2019.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)

#############################
##2020 DATA
Raw2020data <- read.csv(paste(seed_path,
                       "SortedSeeds_raw_MORA_2020.csv",sep=""), header=TRUE)
#rm country column
hdr <- dimnames(Raw2020data)[[2]]
Raw2020data <- Raw2020data[,-(which(hdr=="country"))]

#now change species to codes
Raw2020data$spp[Raw2020data$spp[]=="Abies amabalis"] <- "ABAM"
Raw2020data$spp[Raw2020data$spp[]=="Abies lasiocarpa"] <- "ABLA"
Raw2020data$spp[Raw2020data$spp[]=="Abies procera"] <- "ABPR"
Raw2020data$spp[Raw2020data$spp[]=="Chamaecyparis nootkatensis"] <- "CANO"
Raw2020data$spp[Raw2020data$spp[]=="Pinus contorta"] <- "PICO"
Raw2020data$spp[Raw2020data$spp[]=="Pseudotsuga menziesii"] <- "PSME"
Raw2020data$spp[Raw2020data$spp[]=="Taxus brevifolia"] <- "TABR"
Raw2020data$spp[Raw2020data$spp[]=="Thuja plicata"] <- "THPL"
Raw2020data$spp[Raw2020data$spp[]=="Tsuga heterophylla"] <- "TSHE"
Raw2020data$spp[Raw2020data$spp[]=="Tsuga mertensiana"] <- "TSME"

##now write out
write.csv(Raw2020data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2020.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)

#############################
##2021 DATA
Raw2021data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2021.csv",sep=""), header=TRUE)
#rm country column
hdr <- dimnames(Raw2021data)[[2]]
Raw2021data <- Raw2021data[,-(which(hdr=="country"))]

#first, remove incorrectly spelled T004 and A003
rmrow <- which(Raw2021data$stand=="T004") #remove incorrect spelling, duplicate
Raw2021data <- Raw2021data[-rmrow,]
#first, remove incorrectly spelled T004 and A003
rmrow <- which(Raw2021data$stand=="A003") #remove incorrect spelling, duplicate
Raw2021data <- Raw2021data[-rmrow,]

#now, remove duplicated rows
rmrow <- duplicated(Raw2021data)
Raw2021data <- Raw2021data[rmrow=="FALSE",]

#quick check - still all stand / plot expected?
stnds <- unique(Raw2021data$stand)
for(i in 1:length(stnds)){
  tmp <- Raw2021data[Raw2021data$stand==stnds[i],]
  print(stnds[i])
  print(unique(tmp$trapno))
} #OK

#remove Abies sp. species
rmrow <- which(Raw2021data$spp=="Abies sp.") #remove where seeds not collected
Raw2021data <- Raw2021data[-rmrow,]

#now change species to codes
Raw2021data$spp[Raw2021data$spp[]=="Abies amabalis"] <- "ABAM"
Raw2021data$spp[Raw2021data$spp[]=="Abies lasiocarpa"] <- "ABLA"
Raw2021data$spp[Raw2021data$spp[]=="Abies procera"] <- "ABPR"
Raw2021data$spp[Raw2021data$spp[]=="Alnus rubra"] <- "ALRU"
Raw2021data$spp[Raw2021data$spp[]=="Alnus viridis"] <- "ALVI"
Raw2021data$spp[Raw2021data$spp[]=="Chamaecyparis nootkatensis"] <- "CANO"
Raw2021data$spp[Raw2021data$spp[]=="Picea engelmannii"] <- "PIEN"
Raw2021data$spp[Raw2021data$spp[]=="Pinus contorta"] <- "PICO"
Raw2021data$spp[Raw2021data$spp[]=="Pseudotsuga menziesii"] <- "PSME"
Raw2021data$spp[Raw2021data$spp[]=="Taxus brevifolia"] <- "TABR"
Raw2021data$spp[Raw2021data$spp[]=="Thuja plicata"] <- "THPL"
Raw2021data$spp[Raw2021data$spp[]=="Tsuga heterophylla"] <- "TSHE"
Raw2021data$spp[Raw2021data$spp[]=="Tsuga mertensiana"] <- "TSME"

##now write out
write.csv(Raw2021data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2021.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)

#############################
##2022 DATA
Raw2022data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2022.csv",sep=""), header=TRUE)
#rm country column
hdr <- dimnames(Raw2022data)[[2]]
Raw2022data <- Raw2022data[,-(which(hdr=="country"))]

#Change AM17 to AM16
Raw2022data$stand[Raw2022data$stand[]=="AM17"] <- "AM16"

#remove Mistletoe
rmrow <- which(Raw2022data$spp=="Mistletoe") #remove this species
Raw2022data <- Raw2022data[-rmrow,]

#remove Abies sp.
rmrow <- which(Raw2022data$spp=="Abies sp.") #remove this species
Raw2022data <- Raw2022data[-rmrow,]

#Change Acer to ACCI
Raw2022data$spp[Raw2022data$spp[]=="Acer"] <- "ACCI"

#remove Abies
rmrow <- which(Raw2022data$spp=="Abies") #remove this species
Raw2022data <- Raw2022data[-rmrow,]

#remove NA in species
rmrow <- which(is.na(Raw2022data$spp)==TRUE) #remove NA species- means 0 seeds
Raw2022data <- Raw2022data[-rmrow,]

#now change species to codes
Raw2022data$spp[Raw2022data$spp[]=="Abies amabalis"] <- "ABAM"
Raw2022data$spp[Raw2022data$spp[]=="Abies lasiocarpa"] <- "ABLA"
Raw2022data$spp[Raw2022data$spp[]=="Abies grandis"] <- "ABGR"
Raw2022data$spp[Raw2022data$spp[]=="Abies procera"] <- "ABPR"
Raw2022data$spp[Raw2022data$spp[]=="Acer sp."] <- "ACCI"
Raw2022data$spp[Raw2022data$spp[]=="Alnus rubra"] <- "ALRU"
Raw2022data$spp[Raw2022data$spp[]=="Alnus viridis"] <- "ALVI"
Raw2022data$spp[Raw2022data$spp[]=="Chamaecyparis nootkatensis"] <- "CANO"
Raw2022data$spp[Raw2022data$spp[]=="Picea engelmannii"] <- "PIEN"
Raw2022data$spp[Raw2022data$spp[]=="Pinus contorta"] <- "PICO"
Raw2022data$spp[Raw2022data$spp[]=="Pseudotsuga menziesii"] <- "PSME"
Raw2022data$spp[Raw2022data$spp[]=="Taxus brevifolia"] <- "TABR"
Raw2022data$spp[Raw2022data$spp[]=="Thuja plicata"] <- "THPL"
Raw2022data$spp[Raw2022data$spp[]=="Tsuga heterophylla"] <- "TSHE"
Raw2022data$spp[Raw2022data$spp[]=="Tsuga mertensiana"] <- "TSME"

##now write out
write.csv(Raw2022data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2022.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)

#############################
##2023 DATA
Raw2023data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2023.csv",sep=""), header=TRUE)
#rm country column
hdr <- dimnames(Raw2023data)[[2]]
Raw2023data <- Raw2023data[,-(which(hdr=="country"))]

#now change species to codes
Raw2023data$spp[Raw2023data$spp[]=="Abies amabalis"] <- "ABAM"
Raw2023data$spp[Raw2023data$spp[]=="Abies lasiocarpa"] <- "ABLA"
Raw2023data$spp[Raw2023data$spp[]=="Abies grandis"] <- "ABGR"
Raw2023data$spp[Raw2023data$spp[]=="Abies procera"] <- "ABPR"
Raw2023data$spp[Raw2023data$spp[]=="Acer sp."] <- "ACCI"
Raw2023data$spp[Raw2023data$spp[]=="Alnus rubra"] <- "ALRU"
Raw2023data$spp[Raw2023data$spp[]=="Alnus viridis"] <- "ALVI"
Raw2023data$spp[Raw2023data$spp[]=="Chamaecyparis nootkatensis"] <- "CANO"
Raw2023data$spp[Raw2023data$spp[]=="Picea engelmannii"] <- "PIEN"
Raw2023data$spp[Raw2023data$spp[]=="Pinus contorta"] <- "PICO"
Raw2023data$spp[Raw2023data$spp[]=="Pseudotsuga menziesii"] <- "PSME"
Raw2023data$spp[Raw2023data$spp[]=="Taxus brevifolia"] <- "TABR"
Raw2023data$spp[Raw2023data$spp[]=="Thuja plicata"] <- "THPL"
Raw2023data$spp[Raw2023data$spp[]=="Tsuga heterophylla"] <- "TSHE"
Raw2023data$spp[Raw2023data$spp[]=="Tsuga mertensiana"] <- "TSME"

##now write out
write.csv(Raw2023data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2023.csv", sep=""), 
                                     row.names=FALSE, col.names=TRUE)
#############################
##2024 DATA
Raw2024data <- read.csv(paste(seed_path,
                              "SortedSeeds_raw_MORA_2024.csv",sep=""), header=TRUE)
#rm country column
hdr <- dimnames(Raw2024data)[[2]]
Raw2024data <- Raw2024data[,-(which(hdr=="country"))]

#remove Mistletoe
rmrow <- which(Raw2024data$spp=="Mistletoe") #remove this species
Raw2024data <- Raw2024data[-rmrow,]


#now change species to codes
Raw2024data$spp[Raw2024data$spp[]=="Abies amabalis"] <- "ABAM"
Raw2024data$spp[Raw2024data$spp[]=="Abies amabalis?"] <- "ABAM"
Raw2024data$spp[Raw2024data$spp[]=="Abies lasiocarpa"] <- "ABLA"
Raw2024data$spp[Raw2024data$spp[]=="Abies lasiocarpa?"] <- "ABLA"
Raw2024data$spp[Raw2024data$spp[]=="Abies grandis"] <- "ABGR"
Raw2024data$spp[Raw2024data$spp[]=="Abies procera"] <- "ABPR"
Raw2024data$spp[Raw2024data$spp[]=="Acer"] <- "ACCI"
Raw2024data$spp[Raw2024data$spp[]=="Alnus rubra"] <- "ALRU"
Raw2024data$spp[Raw2024data$spp[]=="Alnus viridis"] <- "ALVI"
Raw2024data$spp[Raw2024data$spp[]=="Chamaecyparis nootkatensis"] <- "CANO"
Raw2024data$spp[Raw2024data$spp[]=="Picea engelmannii"] <- "PIEN"
Raw2024data$spp[Raw2024data$spp[]=="Pinus contorta"] <- "PICO"
Raw2024data$spp[Raw2024data$spp[]=="Pseudotsuga menziesii"] <- "PSME"
Raw2024data$spp[Raw2024data$spp[]=="Taxus brevifolia"] <- "TABR"
Raw2024data$spp[Raw2024data$spp[]=="Thuja plicata"] <- "THPL"
Raw2024data$spp[Raw2024data$spp[]=="Tsuga heterophylla"] <- "TSHE"
Raw2024data$spp[Raw2024data$spp[]=="Tsuga mertensiana"] <- "TSME"

##now write out
write.csv(Raw2024data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2024.csv", sep=""), 
                                    row.names=FALSE, col.names=TRUE)