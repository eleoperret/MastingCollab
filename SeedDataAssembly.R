########################################################################
#### MORA seed data - pull together info needed for modeling        ####
#### Includes sorted seed data, seed trap size data, germinant data ####
#### Adds zeroes, and indicates where data is missing (NAs)         ####
########################################################################

###Read in all clean sorted seed data
seed_files <- list.files(path="./data/rawdata/sortedseeds/clean&notes")
rm_fl <- grep("error", seed_files) #remove error files
seed_files <- seed_files[-rm_fl]