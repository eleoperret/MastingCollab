####################################################################
### script to clean data files from 2009-2024                    ###
### also used to generate 2009-2019 StandTrapInventory Files     ###
### created by Janneke Hille Ris Lambers, 17.10.2025             ###
####################################################################

###Read in 2009 - 2019 seed trap data
## names of data files

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

##Now make corrections to each file, and write out.