###example of script to crunch stats of BC data.####
###Inputs are from the pipeline available at https://github.com/wrana-lab/SPARSEQ_BARCODES

#Required for this example script:

#run27_BC_PairedAnalysis.xlsx
#run27_R1BC_Table.xlsx
#run27_R2BC_Table_Reordered.xlsx
#which are outputted by the SPARSEQ_BARCODES pipeline when running the files for Run 27 BC copy 1 
#available at GEO accession GSE231415
#BC Copy 2 is also available, at GSE246819
#Other BC data sets have been processed exactly like this.

library("ggplot2")
library("dplyr")
library("openxlsx")
library("stringr")


###Run 27 copy 1 - originally a V0 run####
##import row and col bowtie outputs
pythonout<-read.xlsx("run27_BC_PairedAnalysis.xlsx")
head(pythonout)

#exclude the invalid pairs 
pythonout<-subset(pythonout, pythonout$BC_Well != "invalid")

##take python output and make map for each well of which wells caused contamination
pythonout<-read.xlsx("run27_BC_PairedAnalysis.xlsx")
pythonout<-subset(pythonout, pythonout$BC_Well != "invalid")
#due to dupe HEKs we need to merge ID and well 
pythonout$sample_well<-paste(pythonout$sampleID, pythonout$Actual_Well, sep = "_")

pythonout_totals<-as.data.frame(pythonout %>% 
    group_by(sample_well) %>% 
    summarise(countTotal = sum(Count)))
#merge back
pythonout<-merge(pythonout, pythonout_totals, by = "sample_well", all=T)
head(pythonout, 12)
#get percent 
pythonout$percent_of_counts_within_sample<-pythonout$Count / pythonout$countTotal * 100
pythonout$Actual_Well.y<-NULL

python_allwells<-pythonout
python_allwells$rounded_percent<-round(python_allwells$percent_of_counts_within_sample, 2)


python_allwells$row<-gsub("[0-9]", "",python_allwells$BC_Well )
python_allwells$col<-gsub("[A-Z]", "",python_allwells$BC_Well )
python_allwells$ordered_col<-factor(python_allwells$col, levels=c(1:24))
python_allwells$ordered_row<-factor(python_allwells$row, levels=LETTERS[seq( from = 1, to = 16 )])

python_mainwells<-subset(pythonout, pythonout$Actual_Well == pythonout$BC_Well)
python_mainwells$type<-"control"
python_mainwells[grep("^W", python_mainwells$sampleID),]$type<-"sample"

python_mainwells<-subset(python_mainwells, python_mainwells$type == "sample")
mean_correct_percent<-mean(python_mainwells$percent_of_counts_within_sample)

python_allwells_spread<-read.xlsx("run27_BC_PairedAnalysis.xlsx")
python_allwells_spread
python_allwells_spread<-subset(python_allwells_spread, python_allwells_spread$BC_Well != "invalid")

#need total counts for each BC pair, then we will see the percent of that pair that spread around the plate 
python_allwells_spread_totals<-as.data.frame(python_allwells_spread %>% 
    group_by(BC_Well) %>% 
    summarise(countTotal = sum(Count)))
#merge back
python_allwells_spread<-merge(python_allwells_spread, python_allwells_spread_totals, by = "BC_Well", all=T)
head(python_allwells_spread, 12)
#get percent of BC pair found in each diff well
python_allwells_spread$percent_of_counts_within_sample<-python_allwells_spread$Count / python_allwells_spread$countTotal * 100

python_allwells_spread$rounded_percent<-round(python_allwells_spread$percent_of_counts_within_sample, 2)

python_allwells_spread$row<-gsub("[0-9]", "",python_allwells_spread$Actual_Well )
python_allwells_spread$col<-gsub("[A-Z]", "",python_allwells_spread$Actual_Well )
python_allwells_spread$ordered_col<-factor(python_allwells_spread$col, levels=c(1:24))
python_allwells_spread$ordered_row<-factor(python_allwells_spread$row, levels=LETTERS[seq( from = 1, to = 16 )])

python_allwells_spread<-python_allwells_spread[order(python_allwells_spread$sampleID, python_allwells_spread$Actual_Well),]

####next generate table and summary stats of these above tables
python_allwells_correctPair<-subset(python_allwells, python_allwells$Actual_Well == python_allwells$BC_Well)
python_allwells_wrongPair<-subset(python_allwells, python_allwells$Actual_Well != python_allwells$BC_Well)

python_allwells_spread_correctPair<-subset(python_allwells_spread, python_allwells_spread$Actual_Well == python_allwells_spread$BC_Well)
python_allwells_spread_wrongPair<-python_allwells_spread[,1:4]

#resum to get overall values for each BC Pair
ppython_allwells_spread_wrongPair_totals<-as.data.frame(python_allwells_spread_wrongPair %>% 
    group_by(BC_Well) %>% 
    summarise(countTotal = sum(Count)))
#merge back
python_allwells_spread_wrongPair<-merge(python_allwells_spread_wrongPair, ppython_allwells_spread_wrongPair_totals, by = "BC_Well", all=T)
head(python_allwells_spread_wrongPair, 12)
#get percent of BC pair found in each diff well
python_allwells_spread_wrongPair$percent_of_counts_within_sample<-python_allwells_spread_wrongPair$Count / python_allwells_spread_wrongPair$countTotal * 100
head(python_allwells_spread_wrongPair, 60)

#drop correct BCs
python_allwells_spread_wrongPair<-subset(python_allwells_spread_wrongPair, python_allwells_spread_wrongPair$BC_Well != python_allwells_spread_wrongPair$Actual_Well)
head(python_allwells_spread_wrongPair, 60)
#resum percents
ppython_allwells_spread_wrongPair_total_percents<-as.data.frame(python_allwells_spread_wrongPair %>% 
    group_by(BC_Well) %>% 
    summarise(percentTotal = sum(percent_of_counts_within_sample)))

wells<-read.xlsx("run27_R1BC_Table.xlsx")
wells<-wells$well
ppython_allwells_spread_wrongPair_total_percents_trimmed<-subset(ppython_allwells_spread_wrongPair_total_percents, 
                                                                 ppython_allwells_spread_wrongPair_total_percents$BC_Well %in% c(wells))

ppython_allwells_spread_wrongPair_total_percents_trimmed$Experiment<-"Run 27"

list_of_datasets <- list("Correct BCs in Correct Well %" = python_allwells_correctPair, "Correct BCs in Wrong Well %" = ppython_allwells_spread_wrongPair_total_percents_trimmed)
write.xlsx(list_of_datasets, file = "PythonBC_stats_Run27.xlsx")
write.xlsx(python_allwells_spread_wrongPair, file = "python_allwells_spread_wrongPair_Run27.xlsx")

