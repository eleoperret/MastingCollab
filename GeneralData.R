#Objects I use regularly





library(dplyr)

getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting/data")

#Manually added (to do later)
mapRainier<-read.csv("rawdata/Cleaned_mapping_2017 Rainier.csv", sep=";")

#Species per stand if present 
stands_per_species <- mapRainier %>%
  group_by(species) %>%                   # group by species
  summarise(stands = list(unique(stand_id)),  # list of unique stands
            .groups = "drop")             # ungroup

# Save the object
saveRDS(stands_per_species, file = "stands_per_species.rds")


#Stand elevation
stand_elev <- data.frame(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13",
            "AO03","AG05","AV06","AX15","AB08","PP17",
            "AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600,668,700,850,1450,850,
                900,950,1060,1090,1100,1150,
                1150,1200,1450,1600,1700,1800)
)
saveRDS(stand_elev, file= "stand_elevation_table.rds")
