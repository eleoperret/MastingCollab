##############################################################
#### Script to analyze seed masting data from Mt. Rainier ####
#### Script created 12/01/2024; last worked on 15/01/2024 ####
##############################################################

#To do
#Add 2018-2023 data (JHRL working on cleaning and merging scripts)
#Some issues with fitting most complicated model - probably b/c some years with 0 seeds
#Assess how to decide which stands to keep (what is enough data?)
#add basal area (per stand) to analysis (how much do stands vary by year)
#how to identify masting? 

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
#spp <- unique(seeds_anal$species) #all species
spp <- c("ABAM","ABLA","CANO","PICO","PSME","THPL", "TSHE", "TMSE")
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
    }}

    if(sum(na.omit(spp_seeds_anal$filledseeds))<10){next}
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
    AICs <- AIC(null.mod,yr.mod,stnd.mod,yrstnd.mod,yrxstnd.mod)
    print(anova(null.mod,yr.mod,stnd.mod,yrstnd.mod,yrxstnd.mod))
    print(AICs)
    
    #Graph
    # Specify the directory where PNG files will be stored
    output_directory <- "C:/Users/eleop/polybox/phD/PhD/R/Masting_US/Masting/Output/"
    # Create a PNG file for each species
    png(file = paste0(output_directory, "seed_density_", spp[i], ".png"), res = 100)
    #average
    seed_avg <- tapply(spp_seeds_anal$filledseeds, 
                       list(spp_seeds_anal$year,spp_seeds_anal$stand),
                       mean)
    #create a blank plot to add points to
    yrs <- as.numeric(dimnames(seed_avg)[[1]])
    plot(yrs,seed_avg[,1], xlab="years",ylab="seed density",
         type="n", ylim=c(0,1.1*max(na.omit(seed_avg[,]))))
    title(spp[i])
    
    #for loop to add lines for each year
    for(j in 1:dim(seed_avg)[2]){
      points(yrs,seed_avg[,j], type="b", pch=21, bg="yellowgreen")
    }
    # Save the plot as PNG
    dev.off()
}
  
#Second one for species but all sites together
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
    }}
  
  if(sum(na.omit(spp_seeds_anal$filledseeds))<10){next}
  #Graph
  # Specify the directory where PNG files will be stored
  output_directory <- "C:/Users/eleop/polybox/phD/PhD/R/Masting_US/Masting/Output/"
  # Create a PNG file for each species
  png(file = paste0(output_directory, "seed_density_2", spp[i], ".png"), res = 100)
  #average
  seed_avg <- tapply(
    spp_seeds_anal$filledseeds,
    spp_seeds_anal$year,
    mean
  )
  #create a blank plot to add points to
  # create a blank plot to add points to
  yrs <- as.numeric(names(seed_avg))
  plot(yrs, seed_avg, xlab = "years", ylab = "seed density",
       type = "n", ylim = c(0, 1.1 * max(na.omit(seed_avg))))
  title(spp[i])
  
  # add points for each year
  points(yrs, seed_avg, type = "b", pch = 21, bg = "yellowgreen")
  
  # Save the plot as PNG
  dev.off()

}






