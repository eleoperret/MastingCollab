# Masting, Seeds and Seedlings (Mt. Rainier)
Mount Rainier Seed and Seedling (germinants)

This repo contains the data collected from 2008 - 2022(ish). Data from 2009-2017 was examined and cleaned, and you can find a full excel file (all years together) with I think some useful meta data in the two files MORA_cleanseeds_2009-2017.xlsx (also the seed trap size / info). The individual csv files for each year are still in the raw data, and my eventual goal is to add a cleaning and merging script so we can just update the data each year.

The MastingAnalysis.R script is a very basic script I created for Eleonore so she could explore the cleaned 2009-2017 data.

Data collected after 2017 hasn't been integrated, so the data are added as they were entered. This means there could be typos, odd columns, or other oddities in those csv files especially. Please let me know which issues arise so I can fix them...

In this repo is the following data (in raw data folder):
1. Sortedseeds csv files includes the number and species of seeds sorted from seed traps, one file per year. Seeds are distinguished by filled (full embryo) vs. unfilled. Number of cones are also included as counts, as well as the number of seeds dissected from cones. It is important to recognize that the zeroes are not included (for example, if zero TSHE seeds were found in a trap, you have to infer that from the fact that no TSHE seeds were found). Note that the data collected from seed traps in 2012 (and labeled 2012) came from seeds produced and dispersed in fall of 2011. Therefore 2012 labeled seed trap data is the input for 2012 labeled germinant data (this gets confusing).
2. Sortedseeds extra csv file for 2011. In 2010 we decided to collect a few seed traps in fall, to see if seed numbers were different. This file is labeled fall 2011 even though it came from fall 2010, because those numbers should be added to 2011 sorted seed csv.   
3. Germinantsintraps. Sometimes seeds germinate in traps. If that happens, the seed coats would either be counted as unfilled or they possibly would blow away. We count the germinants in traps and remove them before collecting the contents, which means each germinant = 1 viable seed and should be added to filled seed totals. Potentially, this results in an overcount if using unfilled + filled for seed production.
4. Trapinfo. The size of the trap in that year. Pretty much the same for the first 7-8 years.
5. 1styearseedlings. Environmental data and number of first year seedlings (germinants) in 1 m2 quadrats adjacent to each seed trap. For some years, this data will include more than one survey for a stand - this would happen if we came out, thought it was early for a germinant census (snow just melted) and then came back later and redid the survey (sometimes without redoing the environment census). I would take the later date survey, or the one with the greater number of seedlings. Note - 2018 seedling data seem not to have been entered (?!?), but I have the hard copy datasheets (yay!) so this will be done soon.
   

