# Masting, Seeds and Seedlings (Mt. Rainier)
Mount Rainier Seed and seedling data

Updated December 7, 2012

This repo contains the seed data collected from 2009 - 2024. 
- The cleaned data is called SeedData_all.csv. It includes all the seed data (from 2009 -2024), with the zeroes inserted for all species (even if the tree species wasn't found in the stand). It does not include entries for traps that were collected, but not sorted (e.g. lost). The script SeedDataAssembly.R (in the main folder) creates this file from the individual files in the Data folder (see below).

The folder called data has a folder in it called rawdata, and a few old cleaned versions of the data (prior to 2017). The raw data folder includes:
1. Sortedseeds csv files includes the number and species of seeds sorted from seed traps, one file per year with the word raw. As of December 7, 2025, there is now a clean&notes folder that includes the cleaned data (all in the same format, typos removed). The cleaning was done by a script in the main folder (called SeedDataCleaning_2009-2024.R). This cleaning was mostly to homogenize column labels, fix species names, etc. A few notes:
- Zeroes are not included in these raw or cleaned data files. For example, if zero TSHE seeds were found in a trap, then there will not be an entry for TSHE for that trap in the raw data files.
- The data collected from seed traps in 2012 (and labeled 2012) came from seeds produced and dispersed in fall of 2011. Therefore seed trap data from 2012 comes from seeds produced in 2011, and is relevant to the 2012 germinant data.   
3. Germinantsintraps. Sometimes seeds germinate in traps. If that happens, the seed coats would either be counted as unfilled or they possibly would blow away. We count the germinants in traps and remove them before collecting the contents, which means each germinant = 1 viable seed and should be added to filled seed totals.  As of December 7, 2025, there is now a clean&notes folder that includes the cleaned data (all in the same format, typos removed). The R script in in the main folder (called GerminantsinTrapsCleaning_2009-2024.R) was used to clean - mostly homogenizing the format, changing typos, and a few minor changes.
4. Trapinfo. This folder includes two types of data. Trapsizes by year, and something called the standtrap inventory files, which includes information on which traps, of the ones collected, are missing because the sample was lost or the seed trap was destroyed. There are traps labeled 1-15, but the stands differ in when these were sampled. Specifically, not all traps were established in the same year. In the early years, only traps 1-6 were established and sampled - eventually in all 18 stands. After 2020, sampling only occurred in 10 of the 18 stands. Between 2020 and 2025, traps 7-15 were installed in these 10 stands. Traps 1-15 were sampled in one year (the year differs slightly by stand), and then only 7-15 are being sampled from now on. 
5. 1styearseedlings. Environmental data and number of first year seedlings (germinants) in 1 m2 quadrats adjacent to each seed trap. Note- this still needs cleaning / updating as of December 2025

Some additional notes and caveats about the data (from JHRL)
1. Trap size varied over the years
2. Sorted seed data: seeds from cones were counted from 2009-2012, and starting in 2013 (until today) sorted to filled and unfilled.
3. Germinants in traps were counted only starting in 2010. There were some years / stands in the early years where datasheets are missing - likely this was because there were no germinants, but I cannot be 100% sure of that. I have so far assumed that if there is no datasheet there were no germinants in traps. This is for: the following stands / years: All stands collected 2009; PARA in 2010; AM16, AR07, AV02, SPRY, & PARA in 2011, AR07, PARA, SUNR, TO11 in 2012. 
   

