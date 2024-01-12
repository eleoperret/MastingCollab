##############################################################
#### Script to analyze seed masting data from Mt. Rainier ####
#### Script created 12/01/2024; last worked on 12/01/2024 ####
##############################################################

###Read in seed data and seed trap size data
###Note - seed data have been cleaned and compiled from individual years
seeds <- read.csv("./data/MORA_cleanseeds_2009-2017.csv", header=TRUE) #seeds
trapsize <- read.csv("./data/MORA_seedtrapinfo_2009-2017.csv", header=TRUE) #trapsizes

###reconfigure trapsize data
trapsize2 <- trapsize[,1:5]
size <- rep(NA, times=dim(trapsize)[1]); year <- size
trapsize2 <- cbind(trapsize[,1:5], year, size)
finaltrap <- c()

for(i in 2009:2017){
  trapsize3 <- trapsize2
  trapsize3[,6] <- rep(i, dim(trapsize2)[1])
  tmpsize <- trapsize[,i-2002]
  trapsize3[,7] <- tmpsize
  finaltrap <- rbind(finaltrap,trapsize3)
}

#Merge size and seed data
seeds_anal <- merge(finaltrap,seeds,by = c("stand","trap","year"))

#Quick analysis by species
spp <- unique(seeds_anal$species)
stnds <- unique(seeds_anal$stand)

for(i in 1:length(spp)){
  spp_seeds_anal <- seeds_anal[seeds_anal$species==spp[i],]

  #remove stands where the species is never found
  #first identify those stands
  stnds_rm <- c()
  
  for(j in 1:length(stnds)){
    stnd_spp_seeds_anal <- spp_seeds_anal[spp_seeds_anal$stand == stnds[j],]
    if(sum(na.omit(stnd_spp_seeds_anal$filledseeds))==0)
      {stnds_rm <- c(stnds_rm,stnds[j])}
  }
  
  #remove stands with no seeds in any years from dataframe
  if(length(stnds_rm)>15){next}
  
  if(length(stnds_rm)>0){
    for(j in 1:length(stnds_rm)){
      spp_seeds_anal <- spp_seeds_anal[spp_seeds_anal$stand!=stnds_rm[j],]
    }

    if(sum(na.omit(spp_seeds_anal$filledseeds))<10){break}
    #Now analyzee - ask whether year and site and their interaction matter
    #use a glm with a Poisson distribution, and an offset (trapsize)
    null.mod <- glm(filledseeds ~ 1, family = poisson, offset = size, 
                    data = spp_seeds_anal)
    yr.mod <- glm(filledseeds ~ year, family = poisson, offset = size, 
                  data = spp_seeds_anal)
    stnd.mod <- glm(filledseeds ~ stand, family = poisson, offset = size, 
                    data = spp_seeds_anal)
    yrstnd.mod <- glm(filledseeds ~ year + stand, family = poisson, 
                      offset = size, data = spp_seeds_anal)
    yrxstnd.mod <- glm(filledseeds ~ year*stand, family = poisson, 
                       offset = size, data = spp_seeds_anal)
    
    #Anova table
    print(spp[i])
    anova(null.mod,yr.mod,stnd.mod,yrstnd.mod,yrxstnd.mod)
    
    #Graph
    
    
  }
  
  
  
  
}
