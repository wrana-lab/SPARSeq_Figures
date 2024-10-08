####Bar plot of pilot samples for 4A#####

pilot<-read.xlsx("v25data_pilotset_forheatmap.xlsx")

pilotvars<-table(pilot$Variant.Details)
pilotvars<-pilotvars[order(pilotvars, decreasing = T)]
pilotvars<-as.data.frame(pilotvars)
colnames(pilotvars)<-c("variant", "count")
pilotvars
#        variant count
# 1           WT  1142
# 5      B.1.1.7     6
# 15  P1/B.1.351     1

# 2        S691S     9 --> SPBS
# 3        N679K     8 --> SPBS
# 4        A684V     7 --> SPBS
# 6        P681H     6 --> SPBS
# 16 P681H+S691S     1 --> SPBS
# 8        I692I     3 --> SPBS
# 9        P681R     2 --> SPBS
# 10       R682R     2 --> SPBS

# 11       A699C     1 --> RDRP
# 12    A702Stop     1 --> RDRP
# 13       C697S     1 --> RDRP

# 14 E484K+Q677H     1 --> RBM + PBS

# 7        E484K     5 --> SRBM
# 17       S494P     1 --> SRBM
# 18       Y489Y     1 --> SRBM

#so fixed order is
orderedYax<-c("A702Stop", "A699C", "C697S", "E484K+Q677H", "P681H+S691S", "I692I", "S691S", "A684V", "R682R", "P681R", "P681H", "N679K", "S494P", "Y489Y", "E484K", "P1/B.1.351", "B.1.1.7", "WT" )
pilotvars$variant_F<-factor(pilotvars$variant, levels=(orderedYax))

pilotvars$category<-"S-Pbs"
pilotvars[pilotvars$variant%in% c("E484K", "Y489Y", "S494P"),]$category<-"S-Rbm"
pilotvars[pilotvars$variant == "E484K+Q677H",]$category<-"S-Rbm + S-Pbs"
pilotvars[pilotvars$variant%in% c("C697S", "A699C", "A702Stop"),]$category<-"RdRP"
pilotvars[pilotvars$variant == "WT",]$category<-"No Mutation"
pilotvars[pilotvars$variant == "B.1.1.7",]$category<-"VoC"
pilotvars[pilotvars$variant == "P1/B.1.351",]$category<-"VoC"

collist<-c("VoC" = "red", "No Mutation" = "black", "S-Pbs" = "blue", "S-Rbm" = "orange", "S-Rbm + S-Pbs" = "gray", "RdRP" = "yellow" )


#generate 2 copies of the graph, because we want to modify it in illustrator to show the lower values

#full version showing the entire WT bar
pilotvrs_full<-ggplot(pilotvars, aes(x = count, y = variant_F, fill = category)) + geom_col(colour="black") + xlab("Count of Pilot Samples") + 
  ylab("Variant") + ggtitle("Pilot Set Variant Distribution") + scale_fill_manual(values=collist) +  theme_bw() +
  scale_x_continuous(breaks=seq(0, 1200, 100), limits=c(0,1200)) +
    theme(axis.title = element_text(size=10), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text = element_text(size = 8),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
pilotvrs_full

#cut version showing lower numbers
pilotvrs_cut<-ggplot(pilotvars, aes(x = count, y = variant_F, fill = category)) + geom_col(colour="black") + xlab("Count of Pilot Samples") + 
  ylab("Variant") + ggtitle("Pilot Set Variant Distribution") + scale_fill_manual(values=collist) +  theme_bw() +
  #scale_x_continuous(breaks=seq(0, 20, 10), limits=c(0,20)) + 
  coord_cartesian(xlim=c(0,10)) +
    theme(axis.title = element_text(size=10), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text = element_text(size = 8),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
pilotvrs_cut

pdf("pilotset_fig4A.pdf", width = 8, height = 6)
grid.arrange(pilotvrs_full, pilotvrs_cut, ncol=2)
dev.off()


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

pdf("VariantHeatmap1_AllVariantsBinarized_datav25_orderedbyDate.pdf", height = 17, width = 17)
hm_4c
dev.off()




####Variant Curves for panel 4D#####
###get each variant and its + alone on x axis, shown as % of that group

re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v25.xlsx")
#fix dates
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")
re_filtered_maindataset$Date.of.collection.Prcsd<-lubridate::parse_date_time(re_filtered_maindataset$Date.of.collection,"ymd")

# #reorder by date 
re_filtered_maindataset<-re_filtered_maindataset[order(re_filtered_maindataset$Date.of.collection.Prcsd, decreasing = F),]
week(re_filtered_maindataset$Date.of.collection.Prcsd)
# https://lubridate.tidyverse.org/reference/round_date.html#rounding-up-date-objects-1
re_filtered_maindataset$rounded_week<-floor_date(re_filtered_maindataset$Date.of.collection.Prcsd, "week")

#get reordered date
re_filtered_maindataset$dateforgraph<-
  paste(month(re_filtered_maindataset$rounded_week),
        day(re_filtered_maindataset$rounded_week),
        year(re_filtered_maindataset$rounded_week), sep = "-" )

#get order then #force order  
countsperroundedweek_list<-unique(re_filtered_maindataset$dateforgraph)
countsperroundedweek_list_modified<-c(countsperroundedweek_list[1:2], "12-27-2020", countsperroundedweek_list[3], "1-10-2021", countsperroundedweek_list[4:118])
re_filtered_maindataset$dateforgraph_factor<-factor(re_filtered_maindataset$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)

###extract colours matching overall plot 
library(scales)
#hex <- hue_pal()(13)
# #alpha      b/g         delta.    eta         mu    omicron    2.75.    2/3       4/5       other     WT        XBB         XBB.1.5
#c("#F8766D" "#E18A00" "#BE9C00" "#8CAB00" "#24B700" "#00BE70" "#00C1AB" "#00BBDA" "#00ACFC" "#8B93FF" "#D575FE" "#F962DD" "#FF65AC")

###WT subset###
wtonly<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == "WT" | re_filtered_maindataset$Variant_Name_New == "WT+" )
wtonly_weekly_tally<-as.data.frame(wtonly %>% group_by(Variant_Name_New, dateforgraph) %>% tally())
colnames(wtonly_weekly_tally)<-c("Variant_Name_New", "dateforgraph", "weeklysum")
#get total per annotation per week [Variant_Name]
wtonly_perweek<-as.data.frame(wtonly %>% 
   group_by(Variant_Name_New) %>% 
   tally())
#merge weekly sum to dataframe for heatmap
wtonly_weekly<-merge(wtonly_perweek, wtonly_weekly_tally, by ="Variant_Name_New", all=T)
###divide by weekly total to get proportion
wtonly_weekly$weekly_proportion<- wtonly_weekly$weeklysum / wtonly_weekly$n * 100
#force order  - need to keep dates between highest and lowest 
countsperroundedweek_list_modified_wt<-countsperroundedweek_list_modified[countsperroundedweek_list_modified %in% wtonly_weekly$dateforgraph]
wtonly_weekly$dateforgraph_factor<-factor(wtonly_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified_wt)
wtonly_weekly$dateforgraph_date<-lubridate::parse_date_time(wtonly_weekly$dateforgraph,"mdy")
as.Date(wtonly_weekly$dateforgraph_date)

pairgraph_WT2<-ggplot(wtonly_weekly, aes(x=as.Date(dateforgraph_date), y=weekly_proportion,
                group = Variant_Name_New, linetype = Variant_Name_New)) + scale_x_date(date_breaks = "1 month", date_labels = "%b %Y" )+
  scale_linetype_manual(values = c("WT" = 1, "WT+" = 2), name = "Variant") + 
  geom_smooth(se=F, colour = "#DF70F8") + xlab("Collection Date") + ylab("% of Total Cases [within WT+/-]") + theme_bw() + #scale_x_discrete(drop=FALSE) +
  ggtitle('SPARseq WT Data: Dec 2020 - March 2023') + 
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
    
pdf("WT_weekly_data_v25Data.pdf", width = 8, height = 4 )
pairgraph_WT2
dev.off()


###alpha
alphaonly<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == "Alpha" | re_filtered_maindataset$Variant_Name_New == "Alpha+" )
alphaonly_weekly_tally<-as.data.frame(alphaonly %>% group_by(Variant_Name_New, dateforgraph) %>% tally())
colnames(alphaonly_weekly_tally)<-c("Variant_Name_New", "dateforgraph", "weeklysum")
#get total per annotation per week [Variant_Name]
alphaonly_perweek<-as.data.frame(alphaonly %>% 
   group_by(Variant_Name_New) %>% 
   tally())
#merge weekly sum to dataframe for heatmap
alphaonly_weekly<-merge(alphaonly_perweek, alphaonly_weekly_tally, by ="Variant_Name_New", all=T)
###divide by weekly total to get proportion
alphaonly_weekly$weekly_proportion<- alphaonly_weekly$weeklysum / alphaonly_weekly$n * 100
#force order  
alphaonly_weekly$dateforgraph_factor<-factor(alphaonly_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
alphaonly_weekly$dateforgraph_date<-lubridate::parse_date_time(alphaonly_weekly$dateforgraph,"mdy")

pairgraph_alpha2<-ggplot(alphaonly_weekly, aes(x=as.Date(dateforgraph_date), y=weekly_proportion,
                group = Variant_Name_New, linetype = Variant_Name_New)) +  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y" )+
  scale_linetype_manual(values = c("Alpha" = 1, "Alpha+" = 2), name = "Variant") + 
  geom_smooth(se=F, colour = "#F8766D") + xlab("Collection Date") + ylab("% of Total Cases [within Alpha+/-]") + theme_bw() + 
  ggtitle('SPARseq Alpha Data: Dec 2020 - March 2023') +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
    
pdf("alpha_weekly_data_v25Data.pdf", width = 8, height = 4 )
pairgraph_alpha2
dev.off()


###delta
deltaonly<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == "Delta" | re_filtered_maindataset$Variant_Name_New == "Delta+" )
deltaonly_weekly_tally<-as.data.frame(deltaonly %>% group_by(Variant_Name_New, dateforgraph) %>% tally())
colnames(deltaonly_weekly_tally)<-c("Variant_Name_New", "dateforgraph", "weeklysum")
#get total per annotation per week [Variant_Name]
deltaonly_perweek<-as.data.frame(deltaonly %>% 
   group_by(Variant_Name_New) %>% 
   tally())
#merge weekly sum to dataframe for heatmap
deltaonly_weekly<-merge(deltaonly_perweek, deltaonly_weekly_tally, by ="Variant_Name_New", all=T)
###divide by weekly total to get proportion
deltaonly_weekly$weekly_proportion<- deltaonly_weekly$weeklysum / deltaonly_weekly$n * 100
#force order  
deltaonly_weekly$dateforgraph_factor<-factor(deltaonly_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
deltaonly_weekly$dateforgraph_date<-lubridate::parse_date_time(deltaonly_weekly$dateforgraph,"mdy")

pairgraph_delta2<-ggplot(deltaonly_weekly, aes(x=as.Date(dateforgraph_date), y=weekly_proportion,
                group = Variant_Name_New, linetype = Variant_Name_New)) +  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y" )+
  scale_linetype_manual(values = c("Delta" = 1, "Delta+" = 2), name = "Variant") + 
  geom_smooth(se=F, colour = "#B79F00") + xlab("Collection Date") + ylab("% of Total Cases [within Delta+/-]") + theme_bw() + 
  ggtitle('SPARseq Delta Data: Dec 2020 - March 2023') +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
    
pdf("delta_weekly_data_v25Data_separate.pdf", width = 8, height = 4 )
pairgraph_delta2
dev.off()


####omicron
omicrononly<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == "Omicron" | re_filtered_maindataset$Variant_Name_New == "Omicron+" )
omicrononly_weekly_tally<-as.data.frame(omicrononly %>% group_by(Variant_Name_New, dateforgraph) %>% tally())
colnames(omicrononly_weekly_tally)<-c("Variant_Name_New", "dateforgraph", "weeklysum")
#get total per annotation per week [Variant_Name]
omicrononly_perweek<-as.data.frame(omicrononly %>% 
   group_by(Variant_Name_New) %>% 
   tally())
#merge weekly sum to dataframe for heatmap
omicrononly_weekly<-merge(omicrononly_perweek, omicrononly_weekly_tally, by ="Variant_Name_New", all=T)
###divide by weekly total to get proportion
omicrononly_weekly$weekly_proportion<- omicrononly_weekly$weeklysum / omicrononly_weekly$n * 100
#force order  
omicrononly_weekly$dateforgraph_factor<-factor(omicrononly_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
omicrononly_weekly$dateforgraph_date<-lubridate::parse_date_time(omicrononly_weekly$dateforgraph,"mdy")

pairgraph_omicron2<-ggplot(omicrononly_weekly, aes(x=as.Date(dateforgraph_date), y=weekly_proportion,
                group = Variant_Name_New, linetype = Variant_Name_New)) +  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y" )+
  scale_linetype_manual(values = c("Omicron" = 1, "Omicron+" = 2), name = "Variant") + 
  geom_smooth(se=F, colour = "#00BC56") + xlab("Collection Date") + ylab("% of Total Cases [within Omicron+/-]") + theme_bw() +
  ggtitle('SPARseq Omicron Data: Dec 2020 - March 2023') +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
    
pdf("omicron_weekly_data_v25Data.pdf", width = 8, height = 4 )
pairgraph_omicron2
dev.off()


###ba 23
omicrononly<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == "Omicron BA.2/3" | re_filtered_maindataset$Variant_Name_New == "Omicron BA.2/3+" )
omicrononly_weekly_tally<-as.data.frame(omicrononly %>% group_by(Variant_Name_New, dateforgraph) %>% tally())
colnames(omicrononly_weekly_tally)<-c("Variant_Name_New", "dateforgraph", "weeklysum")
#get total per annotation per week [Variant_Name]
omicrononly_perweek<-as.data.frame(omicrononly %>% 
   group_by(Variant_Name_New) %>% 
   tally())
#merge weekly sum to dataframe for heatmap
omicrononly_weekly<-merge(omicrononly_perweek, omicrononly_weekly_tally, by ="Variant_Name_New", all=T)
###divide by weekly total to get proportion
omicrononly_weekly$weekly_proportion<- omicrononly_weekly$weeklysum / omicrononly_weekly$n * 100
#force order  
omicrononly_weekly$dateforgraph_factor<-factor(omicrononly_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
omicrononly_weekly$dateforgraph_date<-lubridate::parse_date_time(omicrononly_weekly$dateforgraph,"mdy")

pairgraph_omicron2<-ggplot(omicrononly_weekly, aes(x=as.Date(dateforgraph_date), y=weekly_proportion,
                group = Variant_Name_New, linetype = Variant_Name_New)) +  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y" )+
  scale_linetype_manual(values = c("Omicron BA.2/3" = 1, "Omicron BA.2/3+" = 2), name = "Variant") + 
  geom_smooth(se=F, colour = "#00BFC4") + xlab("Collection Date") + ylab("% of Total Cases [within Omicron BA.2/3+/-]") + theme_bw() + 
  ggtitle('SPARseq Omicron BA.2/3 Data: Dec 2020 - March 2023') +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
    
pdf("omicronBA23_weekly_data_v25Data.pdf", width = 8, height = 4 )
pairgraph_omicron2
dev.off()


####ba 45
omicrononly<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == "Omicron BA.4/5" | re_filtered_maindataset$Variant_Name_New == "Omicron BA.4/5+" )
omicrononly_weekly_tally<-as.data.frame(omicrononly %>% group_by(Variant_Name_New, dateforgraph) %>% tally())
colnames(omicrononly_weekly_tally)<-c("Variant_Name_New", "dateforgraph", "weeklysum")
#get total per annotation per week [Variant_Name]
omicrononly_perweek<-as.data.frame(omicrononly %>% 
   group_by(Variant_Name_New) %>% 
   tally())
#merge weekly sum to dataframe for heatmap
omicrononly_weekly<-merge(omicrononly_perweek, omicrononly_weekly_tally, by ="Variant_Name_New", all=T)
###divide by weekly total to get proportion
omicrononly_weekly$weekly_proportion<- omicrononly_weekly$weeklysum / omicrononly_weekly$n * 100
#force order  
omicrononly_weekly$dateforgraph_factor<-factor(omicrononly_weekly$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)
omicrononly_weekly$dateforgraph_date<-lubridate::parse_date_time(omicrononly_weekly$dateforgraph,"mdy")

pairgraph_omicron2<-ggplot(omicrononly_weekly, aes(x=as.Date(dateforgraph_date), y=weekly_proportion,
                group = Variant_Name_New, linetype = Variant_Name_New)) +  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y" )+
  scale_linetype_manual(values = c("Omicron BA.4/5" = 1, "Omicron BA.4/5+" = 2), name = "Variant") + 
  geom_smooth(se=F, colour = "#00B6EB") + xlab("Collection Date") + ylab("% of Total Cases [within Omicron BA.4/5+/-]") + theme_bw() + 
  ggtitle('SPARseq Omicron BA.4/5 Data: Dec 2020 - April 2023') +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
    
pdf("omicronBA45_weekly_data_v25Data.pdf", width = 8, height = 4 )
pairgraph_omicron2
dev.off()



