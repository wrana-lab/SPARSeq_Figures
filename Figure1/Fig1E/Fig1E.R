####1e bar graph of runs per week.#####

###Plot runs per week/samples per week as measure of consistency##
re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v26.xlsx")
#fix dates
re_filtered_maindataset$Date.of.collection <- as.Date(re_filtered_maindataset$Date.of.collection, 
                                                      format="%m-%d-%Y", origin = "12-30-1899")
###
weekly_runs<-re_filtered_maindataset
head(weekly_runs)
#column 'date' is the date of processing
table(weekly_runs$Date)
#remove 'Processed " and ","
weekly_runs$Date<-gsub("Processed ", "", weekly_runs$Date)
weekly_runs$Date<-gsub(",", "", weekly_runs$Date)
weekly_runs$Date<-gsub("Sept ", "September ", weekly_runs$Date)
weekly_runs$Date<-gsub("Aug ", "August ", weekly_runs$Date)
#convert to date, same as other dates in analysis
weekly_runs[weekly_runs$Date == "February 25 2022",]
#the feb dates look wrong but are legit, they were processed 2 months late
table(weekly_runs$Date)

weekly_runs$Date.Prcsd<-lubridate::parse_date_time(weekly_runs$Date,"mdy")
# #reorder by date
weekly_runs<-weekly_runs[order(weekly_runs$Date.Prcsd, decreasing = F),]
# https://lubridate.tidyverse.org/reference/round_date.html#rounding-up-date-objects-1
#try rounding down to week
weekly_runs$rounded_week<-floor_date(weekly_runs$Date.Prcsd, "week")
table(weekly_runs$rounded_week)
#with this method every date is shifted to the sunday before. so a sat is shifted 6 days to the previous sunday

#get reordered date
weekly_runs$dateforgraph<-
  paste(month(weekly_runs$rounded_week),
        day(weekly_runs$rounded_week),
        year(weekly_runs$rounded_week), sep = "-" )

#get date list
#count of RUNS per week
weekly_runs$V0run
weekly_runs$V1run
weekly_runs$V2run
#need a col for these run IDs
weekly_runs$RunID<-""
weekly_runs[ !is.na(weekly_runs$V2run), ]$RunID<-paste("V2 ", weekly_runs[!is.na(weekly_runs$V2run),]$V2run)
weekly_runs[ !is.na(weekly_runs$V1run) & weekly_runs$RunID == "", ]$RunID<-
  paste("V1 ", weekly_runs[!is.na(weekly_runs$V1run) & weekly_runs$RunID == "",]$V1run)
weekly_runs[ !is.na(weekly_runs$V0) & weekly_runs$RunID == "", ]$RunID<-
  paste("V0 ", weekly_runs[!is.na(weekly_runs$V0) & weekly_runs$RunID == "",]$V0)
table(weekly_runs$RunID)
unique(weekly_runs$RunID) #365 different runs

weekly_runs_perweek<-as.data.frame(weekly_runs %>%
   group_by(dateforgraph) %>%
   summarize(count_distinct = n_distinct(RunID)))

weekly_runs_perweek<- weekly_runs_perweek %>% arrange(mdy(weekly_runs_perweek$dateforgraph))
#there are 2022s in here. they should be at the end
#because we processed the last Dec samples in Jan 2022

## keep versions separated
weekly_runs$RunType<-""
weekly_runs[grep('V0', weekly_runs$RunID),]$RunType<-"V0"
weekly_runs[grep('V1', weekly_runs$RunID),]$RunType<-"V1"
weekly_runs[grep('V2', weekly_runs$RunID),]$RunType<-"V2"
table(weekly_runs$RunType)
weekly_runs[weekly_runs$RunType == "",]

weekly_runs_perweek_separated<-as.data.frame(weekly_runs %>%
   group_by(dateforgraph, RunType ) %>%
   summarize(count_distinct = n_distinct(RunID)))

weekly_runs_perweek_separated<- weekly_runs_perweek_separated %>% arrange(mdy(weekly_runs_perweek_separated$dateforgraph))

weekly_runs_perweek_separated$dateforgraph_2<-lubridate::parse_date_time(weekly_runs_perweek_separated$dateforgraph,"mdy")

#make bar graph
runsperweek<-ggplot(weekly_runs_perweek_separated, aes(dateforgraph_2, count_distinct, fill = RunType)) + geom_col(colour = "black") +
  ggtitle("SPARSeq runs per week") +  ylab("Number of SPARSeq Runs") + #scale_x_discrete(drop=FALSE) +
  xlab("Week") + theme_bw() + scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 6))
pdf("runs_per_week_v26data.pdf", height = 4, width = 15)
runsperweek
dev.off()
#this was stylized in illustrator
