####################################################################
### script to clean data files from 2009-2024                    ###
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

##now write out
write.csv(Raw2010data, file = paste(seed_path,"clean&notes/",
                                    "SortedSeeds_MORA_2010.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


############
##Germinants in traps data files
##Note - 2009 - 2017 data in one format, 2018 onwards in another