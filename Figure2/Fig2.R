
###final new version of main var graph with + incorporated into normal lines#
re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v25.xlsx")
#fix dates
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")

re_filtered_maindataset$Date.of.collection.Prcsd<-lubridate::parse_date_time(re_filtered_maindataset$Date.of.collection,"ymd")
# #reorder by date 
re_filtered_maindataset<-re_filtered_maindataset[order(re_filtered_maindataset$Date.of.collection.Prcsd, decreasing = F),]
week(re_filtered_maindataset$Date.of.collection.Prcsd)
# https://lubridate.tidyverse.org/reference/round_date.html#rounding-up-date-objects-1
#try rounding down to week
re_filtered_maindataset$rounded_week<-floor_date(re_filtered_maindataset$Date.of.collection.Prcsd, "week")
# #reorder by date 
re_filtered_maindataset<-re_filtered_maindataset[order(re_filtered_maindataset$Date.of.collection.Prcsd, decreasing = F),]
# https://lubridate.tidyverse.org/reference/round_date.html#rounding-up-date-objects-1
#try rounding down to week
re_filtered_maindataset$rounded_week<-floor_date(re_filtered_maindataset$Date.of.collection.Prcsd, "week")
#get reordered date
re_filtered_maindataset$dateforgraph<-
  paste(month(re_filtered_maindataset$rounded_week),
        day(re_filtered_maindataset$rounded_week),
        year(re_filtered_maindataset$rounded_week), sep = "-" )
re_filtered_maindataset$dateforgraph
#get order then #force order  
countsperroundedweek_list<-unique(re_filtered_maindataset$dateforgraph)
countsperroundedweek_list_modified<-c(countsperroundedweek_list[1:2], "12-27-2020", countsperroundedweek_list[3], "1-10-2021", countsperroundedweek_list[4:118])
re_filtered_maindataset$dateforgraph_factor<-factor(re_filtered_maindataset$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
tail(re_filtered_maindataset,10)

#get counts per week/proportion per week
#this is the weekly count of SAMPLES 
re_filtered_maindataset_weeklytally<-as.data.frame(re_filtered_maindataset %>% group_by(dateforgraph) %>% tally())
colnames(re_filtered_maindataset_weeklytally)<-c("dateforgraph", "weeklysum")

#get total per annotation per week - this time all + need to be one
re_filtered_maindatasetplustogether<-re_filtered_maindataset
re_filtered_maindatasetplustogether$Variant_Name_New_Trimmed<-gsub("\\+", "", re_filtered_maindatasetplustogether$Variant_Name_New)
  
re_filtered_maindatasetplustogether_perweek<-as.data.frame(re_filtered_maindatasetplustogether %>% 
   group_by(Variant_Name_New_Trimmed, dateforgraph) %>% 
   tally())

#merge weekly sum to dataframe for heatmap
re_filtered_maindatasetplustog_weekly<-merge(re_filtered_maindatasetplustogether_perweek, re_filtered_maindataset_weeklytally, by ="dateforgraph", all=T)
head(re_filtered_maindatasetplustog_weekly)

###divide by weekly total to get proportion
re_filtered_maindatasetplustog_weekly$weekly_proportion<- re_filtered_maindatasetplustog_weekly$n / re_filtered_maindatasetplustog_weekly$weeklysum * 100
head(re_filtered_maindatasetplustog_weekly)

#force order  
re_filtered_maindatasetplustog_weekly$dateforgraph_factor<-factor(re_filtered_maindatasetplustog_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
tail(re_filtered_maindatasetplustog_weekly,10)
head(re_filtered_maindatasetplustog_weekly)

#put together 
overallplot1Av3<-ggplot(re_filtered_maindatasetplustog_weekly, aes(x=dateforgraph_factor, y=n, colour = Variant_Name_New_Trimmed, group = Variant_Name_New_Trimmed)) +
  geom_line(size = .75) + ylab("Weekly Case Count") + theme_bw() + labs(colour = "Variant") + scale_x_discrete(drop=FALSE) +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(hjust = 1, size = 8))

pdf("combined_weekly_data_v25Data_PlusIncorporatedIntoNormalLines.pdf", width = 18, height = 6 )
overallplot1Av3
dev.off()

#####make a main var graph which caps at 1100 to show the patterns without omicron spike to make the figure clearer
max(re_filtered_maindatasetplustog_weekly$n) #6475
re_filtered_maindatasetplustog_weekly[re_filtered_maindatasetplustog_weekly$n > 1200,]
re_filtered_maindatasetplustog_weekly[re_filtered_maindatasetplustog_weekly$dateforgraph == "1-2-2022",]

overallplot1_capped2k<-ggplot(re_filtered_maindatasetplustog_weekly, aes(x=dateforgraph_factor, y=n,
      colour = Variant_Name_New_Trimmed, group = Variant_Name_New_Trimmed)) +  geom_line(size = .75) + ylab("Weekly Case Count [capped at 1100]") + theme_bw() +
  labs(colour = "Variant") + scale_x_discrete(drop=FALSE) + coord_cartesian(ylim=c(0,1000)) + scale_y_continuous(breaks = seq(0,1000, 100)) +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1), axis.title.x=element_blank(),
      axis.text.y = element_text(hjust = 1, size = 8))
overallplot1_capped2k

pdf("combined_weekly_data_v25_Cappedat1000.pdf", width = 18, height = 6 )
overallplot1_capped2k
dev.off()
