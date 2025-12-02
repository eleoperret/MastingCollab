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
##

