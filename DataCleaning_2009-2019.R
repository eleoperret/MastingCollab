####################################################################
### script to clean data files from 2009-2019                    ###
### also used to generate 2009-2019 StandTrapInventory Files     ###
### created by Janneke Hille Ris Lambers, 17.10.2025             ###
####################################################################

###Read in 2009 - 2019 seed trap data
## names of data files

seed_files <- list.files(path="./data/rawdata/sortedseeds/2009-2019_uncleaned")
seed_path <- "./data/rawdata/sortedseeds/2009-2019_uncleaned/"
outputfile <- paste(seed_path,"errorcheck_10.2025.txt", sep="")
write("October 17, 2025 ERROR CHECK", outputfile)

for(i in 2:length(seed_files)){
  tempseed <- read.csv(paste(seed_path, seed_files[i], sep=""))
  #Signal new file
  write("NEW FILE CHECK", file=outputfile, append=TRUE)
  write(seed_files[i], file=outputfile, append=TRUE)
  
  #check years
  write("years", file=outputfile, append=TRUE)
  write(unique(tempseed$Year), file=outputfile, append=TRUE)
  #check stand
  write("stands", file=outputfile, append=TRUE)
  write(unique(tempseed$Stand), file=outputfile, append=TRUE)
  #check species
  write("species", file=outputfile, append=TRUE)
  write(unique(tempseed$Species), file=outputfile, append=TRUE)
  
  ###Now check stand / traps present in that year
  Stands <- unique(tempseed$Stand)
  write("Check traps per stand", file=outputfile, append=TRUE)
  
  for(j in 1:length(Stands)){
    tempseed_stand <- tempseed[tempseed$Stand==Stands[j],]
    write(Stands[j], file=outputfile, append=TRUE)
    write("traps represent", file=outputfile, append=TRUE)
    write(unique(tempseed_stand$Trapno), file=outputfile, append=TRUE)
  }
}

##Now make corrections to each file, and write out.