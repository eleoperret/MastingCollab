####################################################################
### script to clean germinants in trap files from 2009-2024      ###
### also used to generate 2009-2019 StandTrapInventory Files     ###
### created by Janneke Hille Ris Lambers, 02.12.2025             ###
####################################################################

###Read in 2009 - 2024 germinant data
##Note - 2009-2017 one format, 2018 - 2024 different format
###Check species, NAs, etc
###creates errorcheck file
###separate text file provides info on solutions

##Germinant data
germ_files <- list.files(path="./data/rawdata/germinantspertrap")
germ_files <- germ_files[-1]
germ_path <- "./data/rawdata/germinantspertrap/"
outputfile <- paste(germ_path,"clean&notes/errorcheck_12.2025.txt", sep="")
write("December 2, 2025 ERROR CHECK", outputfile)
cat("\n", file = outputfile, append = TRUE) # Adds the blank line

for(i in 1:length(germ_files)){
  tempgerm <- read.csv(paste(germ_path, germ_files[i], sep=""))
  
  #extract year
  tmp1 <- strsplit(germ_files[i],".csv")[[1]]
  yr <- strsplit(tmp1,"Germinants_")[[1]][2]
  yr <- as.numeric(yr)
  
  #Signal new file
  write("NEW FILE CHECK", file=outputfile, append=TRUE)
  write(germ_files[i], file=outputfile, append=TRUE)
  cat("\n", file = outputfile, append = TRUE) # Adds the blank line
  
  #check years
  write("years", file=outputfile, append=TRUE)
  write(unique(tempgerm$Year), file=outputfile, append=TRUE)
  cat("\n", file = outputfile, append = TRUE) # Adds the blank line
  
  #check stand
  write("stands", file=outputfile, append=TRUE)
  write(unique(tempgerm$Stand), file=outputfile, append=TRUE)
  cat("\n", file = outputfile, append = TRUE) # Adds the blank line
  
  #check species, if years 2018 onward
  if(yr>2017){
    write("species", file=outputfile, append=TRUE)
    write(unique(tempgerm$Species), file=outputfile, append=TRUE)
    cat("\n", file = outputfile, append = TRUE) # Adds the blank line
  }
  
  ###Now check stand / traps present in that year
  Stands <- unique(tempgerm$Stand)
  write("Check traps per stand", file=outputfile, append=TRUE)
  
  for(j in 1:length(Stands)){
    tempgerm_stand <- tempgerm[tempgerm$Stand==Stands[j],]
    write(Stands[j], file=outputfile, append=TRUE)
    write("traps present", file=outputfile, append=TRUE)
    write(unique(tempgerm_stand$Trapno), file=outputfile, append=TRUE)
    
    #What about germinant entries
    if(yr<2018){allgerms <- c(unlist(tempgerm_stand[,4:10]))}
    if(yr>2017){allgerms <- tempgerm_stand$Germinants}   
    write("germinant entries", file=outputfile, append=TRUE)  
    write(unique(allgerms), file=outputfile, append=TRUE) 
    cat("\n", file = outputfile, append = TRUE) # Adds the blank line
  }
  cat("\n", file = outputfile, append = TRUE) # Adds the blank line
  cat("\n", file = outputfile, append = TRUE) # Adds the blank line
}

##Now make corrections to individual germinants data, write to clean folder
##Now make corrections to individual sorted seed data, write to clean folder
##

#############################
#############################
##2009 DATA
Raw2009data <- read.csv(paste(germ_path,
                              "Germinants_2009.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2009data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2009data)[[2]] <- hdr2

Clean2009data <- c()

for(i in 4:10){
  tmpdata <- Raw2009data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2009data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2009data <- rbind(Clean2009data,tmpdata2)
}

#remove NAs
Clean2009data2 <- Clean2009data[na.omit(Clean2009data[,5]),]

#Write out file
write.csv(Clean2009data2, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2009.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)

#############################
#############################
##2010 DATA
Raw2010data <- read.csv(paste(germ_path,
                              "Germinants_2010.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2010data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2010data)[[2]] <- hdr2

Clean2010data <- c()

for(i in 4:10){
  tmpdata <- Raw2010data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2010data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2010data <- rbind(Clean2010data,tmpdata2)
}

#remove NAs
Clean2010data2 <- Clean2010data[is.na(Clean2010data[,5])==FALSE,]

#remove 0 entries
Clean2010data3 <- Clean2010data2[Clean2010data2[,5]>0,]

#Write out file
write.csv(Clean2010data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2010.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)

#############################
#############################
##2011 DATA
Raw2011data <- read.csv(paste(germ_path,
                              "Germinants_2011.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2011data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2011data)[[2]] <- hdr2

Clean2011data <- c()

for(i in 4:10){
  tmpdata <- Raw2011data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2011data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2011data <- rbind(Clean2011data,tmpdata2)
}

#remove NAs
Clean2011data2 <- Clean2011data[is.na(Clean2011data[,5])==FALSE,]

#remove 0 entries
Clean2011data3 <- Clean2011data2[Clean2011data2[,5]>0,]

#Write out file
write.csv(Clean2011data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2011.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2012 DATA
Raw2012data <- read.csv(paste(germ_path,
                              "Germinants_2012.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2012data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2012data)[[2]] <- hdr2

Clean2012data <- c()

for(i in 4:10){
  tmpdata <- Raw2012data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2012data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2012data <- rbind(Clean2012data,tmpdata2)
}

#remove NAs
Clean2012data2 <- Clean2012data[is.na(Clean2012data[,5])==FALSE,]

#remove 0 entries
Clean2012data3 <- Clean2012data2[Clean2012data2[,5]>0,]

#Write out file
write.csv(Clean2012data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2012.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2013 DATA
Raw2013data <- read.csv(paste(germ_path,
                              "Germinants_2013.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2013data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2013data)[[2]] <- hdr2

Clean2013data <- c()

for(i in 4:10){
  tmpdata <- Raw2013data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2013data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2013data <- rbind(Clean2013data,tmpdata2)
}

#remove NAs
Clean2013data2 <- Clean2013data[is.na(Clean2013data[,5])==FALSE,]

#remove 0 entries
Clean2013data3 <- Clean2013data2[Clean2013data2[,5]>0,]

#Write out file
write.csv(Clean2013data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2013.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2014 DATA
Raw2014data <- read.csv(paste(germ_path,
                              "Germinants_2014.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2014data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2014data)[[2]] <- hdr2

Clean2014data <- c()

for(i in 4:10){
  tmpdata <- Raw2014data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2014data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2014data <- rbind(Clean2014data,tmpdata2)
}

#remove NAs
Clean2014data2 <- Clean2014data[is.na(Clean2014data[,5])==FALSE,]

#remove 0 entries
Clean2014data3 <- Clean2014data2[Clean2014data2[,5]>0,]

#Write out file
write.csv(Clean2014data3, file = paste(germ_path,"clean&notes/",
                                       "Germinants_MORA_2014.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2015 DATA
Raw2015data <- read.csv(paste(germ_path,
                              "Germinants_2015.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2015data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2015data)[[2]] <- hdr2

Clean2015data <- c()

for(i in 4:10){
  tmpdata <- Raw2015data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2015data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2015data <- rbind(Clean2015data,tmpdata2)
}

#remove NAs
Clean2015data2 <- Clean2015data[is.na(Clean2015data[,5])==FALSE,]

#remove 0 entries
Clean2015data3 <- Clean2015data2[Clean2015data2[,5]>0,]

#Write out file
write.csv(Clean2015data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2015.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2016 DATA
Raw2016data <- read.csv(paste(germ_path,
                              "Germinants_2016.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2016data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2016data)[[2]] <- hdr2

Clean2016data <- c()

for(i in 4:10){
  tmpdata <- Raw2016data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2016data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2016data <- rbind(Clean2016data,tmpdata2)
}

#remove NAs
Clean2016data2 <- Clean2016data[is.na(Clean2016data[,5])==FALSE,]

#remove 0 entries
Clean2016data3 <- Clean2016data2[Clean2016data2[,5]>0,]

#Write out file
write.csv(Clean2016data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2016.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2017 DATA
Raw2017data <- read.csv(paste(germ_path,
                              "Germinants_2017.csv",sep=""), header=TRUE)
#change headers
hdr <- dimnames(Raw2017data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
dimnames(Raw2017data)[[2]] <- hdr2

Clean2017data <- c()

for(i in 4:10){
  tmpdata <- Raw2017data[,c(1:3,i)]
  tmpspp <- rep(dimnames(Raw2017data)[[2]][i], times=dim(tmpdata)[1])
  tmpdata2 <- cbind(tmpdata[,1:3],tmpspp, tmpdata[,4])
  dimnames(tmpdata2)[[2]][4] <- "spp"  
  dimnames(tmpdata2)[[2]][5] <- "germintrap"
  Clean2017data <- rbind(Clean2017data,tmpdata2)
}

#remove NAs
Clean2017data2 <- Clean2017data[is.na(Clean2017data[,5])==FALSE,]

#remove 0 entries
Clean2017data3 <- Clean2017data2[Clean2017data2[,5]>0,]

#Write out file
write.csv(Clean2017data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2017.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)



#############################
#############################
##2018 DATA
Raw2018data <- read.csv(paste(germ_path,
                              "Germinants_2018.csv",sep=""), header=TRUE)

#remove extra rows
rwsin <- which(Raw2018data$Year==2018)
Raw2018data <- Raw2018data[rwsin,]

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2018data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2018data)[[2]]=="Field.Notes")
rm3 <- which(dimnames(Raw2018data)[[2]]=="Data.Entry.Notes")

Raw2018data <- Raw2018data[,-c(rm1,rm2,rm3)]

#change headers
hdr <- dimnames(Raw2018data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2018data)[[2]] <- hdr2

#remove NAs
Clean2018data <- Raw2018data[is.na(Raw2018data$spp)==FALSE,]

#remove blank entry
Clean2018data2 <- Clean2018data[Clean2018data$spp!="",]

#remove 0 entries
Clean2018data3 <- Clean2018data2[Clean2018data2[,5]>0,]

#Write out file
write.csv(Clean2018data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2018.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2019 DATA
Raw2019data <- read.csv(paste(germ_path,
                              "Germinants_2019.csv",sep=""), header=TRUE)

#remove extra rows
rwsin <- which(Raw2019data$Year==2019)
Raw2019data <- Raw2019data[rwsin,]

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2019data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2019data)[[2]]=="Field.Notes")
rm3 <- which(dimnames(Raw2019data)[[2]]=="Data.Entry.Notes")

Raw2019data <- Raw2019data[,-c(rm1,rm2,rm3)]

#change headers
hdr <- dimnames(Raw2019data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2019data)[[2]] <- hdr2

#remove NAs
Clean2019data <- Raw2019data[is.na(Raw2019data$spp)==FALSE,]

#remove blank entry, unknown entry
Clean2019data2 <- Clean2019data[Clean2019data$spp!="-",]
Clean2019data2 <- Clean2019data2[Clean2019data2$spp!="unknown",]

#remove zero entries if any
Clean2019data3 <- Clean2019data2[Clean2019data2[,5]>0,]

#Write out file
write.csv(Clean2019data3, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2019.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2020 DATA
Raw2020data <- read.csv(paste(germ_path,
                              "Germinants_2020.csv",sep=""), header=TRUE)

#remove extra rows
rwsin <- which(Raw2020data$Year==2020)
Raw2020data <- Raw2020data[rwsin,]

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2020data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2020data)[[2]]=="Field.Notes")
rm3 <- which(dimnames(Raw2020data)[[2]]=="Data.Entry.Notes")

Raw2020data <- Raw2020data[,-c(rm1,rm2,rm3)]

#change headers
hdr <- dimnames(Raw2020data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2020data)[[2]] <- hdr2

#remove NAs
Clean2020data <- Raw2020data[is.na(Raw2020data$spp)==FALSE,]

#Write out file
write.csv(Clean2020data, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2020.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2021 DATA
Raw2021data <- read.csv(paste(germ_path,
                              "Germinants_2021.csv",sep=""), header=TRUE)

#remove extra rows
rwsin <- which(Raw2021data$Year==2021)
Raw2021data <- Raw2021data[rwsin,]

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2021data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2021data)[[2]]=="Field.Notes")
rm3 <- which(dimnames(Raw2021data)[[2]]=="Data.Entry.Notes")

Raw2021data <- Raw2021data[,-c(rm1,rm2,rm3)]

#change headers
hdr <- dimnames(Raw2021data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2021data)[[2]] <- hdr2

#remove spp=none
Clean2021data <- Raw2021data[Raw2021data$spp!="none",]

#Write out file
write.csv(Clean2021data, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2021.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)

#############################
#############################
##2022 DATA
Raw2022data <- read.csv(paste(germ_path,
                              "Germinants_2022.csv",sep=""), header=TRUE)

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2022data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2022data)[[2]]=="Field.Notes")
rm3 <- which(dimnames(Raw2022data)[[2]]=="Data.Entry.Notes")

Raw2022data <- Raw2022data[,-c(rm1,rm2,rm3)]

#change headers
hdr <- dimnames(Raw2022data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2022data)[[2]] <- hdr2

#remove spp=none
Clean2022data <- Raw2022data[Raw2022data$spp!="",]

#replace Ab. sp. with ABAM (AE10)
Clean2022data$spp[Clean2022data$spp=="Ab. sp."] <- "ABAM"

# remove rows with AB SP.
Clean2022data2 <- Clean2022data[Clean2022data$spp!="AB SP.",]

#Write out file
write.csv(Clean2022data2, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2022.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2023 DATA
Raw2023data <- read.csv(paste(germ_path,
                              "Germinants_2023.csv",sep=""), header=TRUE)

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2023data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2023data)[[2]]=="Field.Notes")
rm3 <- which(dimnames(Raw2023data)[[2]]=="Data.Entry.Notes")

Raw2023data <- Raw2023data[,-c(rm1,rm2,rm3)]

#change headers
hdr <- dimnames(Raw2023data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2023data)[[2]] <- hdr2

#remove NAs
Clean2023data <- Raw2023data[is.na(Raw2023data$spp)==FALSE,]

#replace Tssp with TSHE (note indicates its likely this species)
Clean2023data$spp[Clean2023data$spp=="TSSP"] <- "TSHE"

#Write out file
write.csv(Clean2023data, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2023.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)


#############################
#############################
##2024 DATA
Raw2024data <- read.csv(paste(germ_path,
                              "Germinants_2024.csv",sep=""), header=TRUE)

##Remove Date, Field.Notes, Data.Entry.Notes columns
rm1 <- which(dimnames(Raw2024data)[[2]]=="Date")
rm2 <- which(dimnames(Raw2024data)[[2]]=="Field.Notes")

Raw2024data <- Raw2024data[,-c(rm1,rm2)]

#change headers
hdr <- dimnames(Raw2024data)[[2]]
hdr2 <- hdr
hdr2[hdr=="Year"]<-"year"
hdr2[hdr=="Stand"]<-"stand"
hdr2[hdr=="Trapno"]<-"trapno"
hdr2[hdr=="Species"]<-"spp"
hdr2[hdr=="Germinants"]<-"germintrap"
dimnames(Raw2024data)[[2]] <- hdr2

#remove NAs
Clean2024data <- Raw2024data[Raw2024data$germintrap>0,]


#Write out file
write.csv(Clean2024data, file = paste(germ_path,"clean&notes/",
          "Germinants_MORA_2024.csv", sep=""), 
          row.names=FALSE, col.names=TRUE)





