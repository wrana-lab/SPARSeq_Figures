
###Binarized Heatmap For Fig 4C####

re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v25.xlsx")
#fix dates
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")

data_for_hm<-re_filtered_maindataset

data_for_hm$Date.of.collection.Prcsd<-lubridate::parse_date_time(data_for_hm$Date.of.collection,"ymd")
# #reorder by date 
data_for_hm<-data_for_hm[order(data_for_hm$Date.of.collection.Prcsd, decreasing = F),]
week(data_for_hm$Date.of.collection.Prcsd)
# https://lubridate.tidyverse.org/reference/round_date.html#rounding-up-date-objects-1
#try rounding down to week
data_for_hm$rounded_week<-floor_date(data_for_hm$Date.of.collection.Prcsd, "week")
#with this method every date is shifted to the sunday before.
#so a sat is shifted 6 days to the previous sunday

# rounding down to week
data_for_hm$rounded_week<-floor_date(data_for_hm$Date.of.collection.Prcsd, "week")
#get reordered date
data_for_hm$dateforgraph<-
  paste(month(data_for_hm$rounded_week),
        day(data_for_hm$rounded_week),
        year(data_for_hm$rounded_week), sep = "-" )
data_for_hm$dateforgraph
#get order then #force order  
countsperroundedweek_list<-unique(data_for_hm$dateforgraph)
countsperroundedweek_list_modified<-c(countsperroundedweek_list[1:2], "12-27-2020", countsperroundedweek_list[3], "1-10-2021", countsperroundedweek_list[4:118])
#counts per week
data_for_hm$Variant.Details
data_for_hm_perweek<-as.data.frame(data_for_hm %>% 
   group_by(dateforgraph) %>% 
   dplyr::count(Variant.Details))

#force order  
data_for_hm_perweek$dateforgraph_factor<-factor(data_for_hm_perweek$dateforgraph, ordered=T, levels=countsperroundedweek_list)

#get date list
#get proportion per week
data_for_hm_perweek_weeklysum<-as.data.frame(data_for_hm %>% group_by(dateforgraph) %>% dplyr::count())
colnames(data_for_hm_perweek_weeklysum)<-c("dateforgraph", "weeklysum")
#merge weekly sum to dataframe for heatmap
data_for_hm_perweek<-merge(data_for_hm_perweek, data_for_hm_perweek_weeklysum, by ="dateforgraph", all=T)

#do with binary version
data_for_hm_perweek_binary<-data_for_hm_perweek
data_for_hm_perweek_binary$binary<-0
data_for_hm_perweek_binary[data_for_hm_perweek_binary$n>0,]$binary<-1
table(data_for_hm_perweek_binary$binary) #all 1

#redo binary version but with rows ordered by first date of occurrence
data_for_hm_perweek_binary_ordered<-data_for_hm_perweek_binary

data_for_hm_perweek_binary_ordered_slim<-data_for_hm_perweek_binary_ordered[,c("dateforgraph","Variant.Details") ]
data_for_hm_perweek_binary_ordered_slim$dateforgraph <- lubridate::mdy(data_for_hm_perweek_binary_ordered_slim$dateforgraph)

#start dates
mindatelist <- as.data.frame(data_for_hm_perweek_binary_ordered_slim %>% 
  dplyr::group_by(Variant.Details) %>% 
  dplyr::mutate(
    date_min = min(dateforgraph)  )) 

unique(mindatelist$Variant.Details) 
#drop dateforgraph and take unique
mindatelist$dateforgraph<-NULL
mindatelist<-unique(mindatelist)
nrow(mindatelist)
#order properly
mindatelist$date_min <- lubridate::ymd(mindatelist$date_min)
mindatelist<-dplyr::arrange(mindatelist, date_min)    

orderofdates<-mindatelist$Variant.Details
revorderofdates<-rev(orderofdates)
#apply order
data_for_hm_perweek_binary_ordered$ordered_date<-factor(data_for_hm_perweek_binary_ordered$Variant.Details, levels = revorderofdates)

hm_4c<-ggplot(data_for_hm_perweek_binary_ordered, aes(x=dateforgraph_factor, y=ordered_date, fill = (binary) )) + 
  ggtitle("Dec 2020 to March 2023: All Variants, Binarized [ordered by date of first appearance]") +
  geom_tile() + xlab("Week") + ylab("Variant") + labs(fill="Presence of Variant") + theme_bw() +
  scale_fill_gradientn(colors=c(low = 'white',  high = 'black'), limits=c(0,1), breaks=c(0,1) )+
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), 
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 4)  ) 
hm_4c

hm_4c<-ggplot(data_for_hm_perweek_binary_ordered, aes(x=dateforgraph_factor, y=ordered_date, fill = (binary) )) + 
  ggtitle("Dec 2020 to March 2023: All Variants, Binarized [ordered by date of first appearance]") +
  geom_tile() + xlab("Week") + ylab("Variant") + labs(fill="Presence of Variant") + theme_bw() +
  scale_fill_gradientn(colors=c(low = 'white',  high = 'black'), limits=c(0,1), breaks=c(0,1) )+
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), 
      axis.text.x = element_blank(), axis.text.y = element_blank())
hm_4c

pdf("VariantHeatmap1_AllVariantsBinarized_datav25_orderedbyDate.pdf", height = 17, width = 17)
hm_4c
hm_4c
dev.off()
