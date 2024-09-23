
###Script to produce Fig 3 analysis and plots from SPARseq manuscript.####
##Revisions version, Sep 2024. 

###load libraries
options(scipen = 9999)
library(openxlsx)
library(lubridate)
library(dplyr)
library(ggplot2)

##This version of analysis uses proportion of sample per day to eliminate effects of different number of samples per day#
#start from main dataset:

re_filtered_maindataset<-read.xlsx("../filtered_data_for_paper_v23_June2024_CleanedVariants_GithubVersion_LC.xlsx")
#adjust dates
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")
re_filtered_maindataset$Variant_Name

data_for_hm<-re_filtered_maindataset
data_for_hm$Date.of.collection.Prcsd<-lubridate::parse_date_time(data_for_hm$Date.of.collection,"ymd")
#reorder by date 
data_for_hm<-data_for_hm[order(data_for_hm$Date.of.collection.Prcsd, decreasing = F),]

#keep using rounded_week
#counts per week
data_for_hm$Variant_Name_New #formerly was VariantDetails
data_for_hm$Variant_Name_NewNoPlus<-data_for_hm$Variant_Name_New 
data_for_hm$Variant_Name_NewNoPlus<-gsub("\\+","",data_for_hm$Variant_Name_NewNoPlus)
data_for_hm_perday<-as.data.frame(data_for_hm %>% 
   group_by(Date.of.collection.Prcsd) %>% 
   dplyr::count(Variant_Name_NewNoPlus))
table(data_for_hm_perday$Variant_Name_NewNoPlus)
head(data_for_hm_perday)
colnames(data_for_hm_perday)<-c("Date.of.collection.Prcsd", "Variant_Name_NewNoPlus",  "VariantSum")

data_for_hm_perdayTotalSum<-as.data.frame(data_for_hm %>% 
   group_by(Date.of.collection.Prcsd) %>% 
   dplyr::count())
head(data_for_hm_perdayTotalSum)
colnames(data_for_hm_perdayTotalSum)<-c("Date.of.collection.Prcsd", "TotalDailySum")

##merge totalsum back to perday and divide to get percent per day 
data_for_hm_perday<-merge(data_for_hm_perday, data_for_hm_perdayTotalSum, by = "Date.of.collection.Prcsd" )
head(data_for_hm_perday, 120)
#basically by early to mid march 2021 Alpha is starting to get over 50% of the daily cases 

data_for_hm_perday$DailyPercent<-data_for_hm_perday$VariantSum/data_for_hm_perday$TotalDailySum * 100

#save file
write.xlsx(data_for_hm_perday, file="data_for_hm_perDAY_PROPORTION_crunchedDataForCompetitions_v24dataset.xlsx")


####organize data####
###based on assessing the file data_for_hm_perDAY_PROPORTION_crunchedDataForCompetitions_v24dataset.xlsx, the transition dates are
#March 5-9 2021 is where WT drops below 50% for 5 straight days
##June 15-20 2021 is where Alpha drops below 50% for 5 straight days. technically one day after that it pops up but w/e
## Dec 13-19 2021 is where Delta drops below 50% 
##We end Omicron at Dec 31 2021 anyway due to changes in government testing policies

#reload and get individual files for each variant

data_for_hm_perday<-read.xlsx("data_for_hm_perDAY_PROPORTION_crunchedDataForCompetitions_v24dataset.xlsx")
data_for_hm_perday$day<-as.Date(data_for_hm_perday$Date.of.collection.Prcsd, origin = "1899-12-30")
#do this as proportion not percent
data_for_hm_perday$DailyPercent<-data_for_hm_perday$DailyPercent/100

cutdata_alpha<-subset(data_for_hm_perday, data_for_hm_perday$Variant_Name_NewNoPlus == "Alpha")
cutdata_wt<-subset(data_for_hm_perday, data_for_hm_perday$Variant_Name_NewNoPlus == "WT")
cutdata_delta<-subset(data_for_hm_perday, data_for_hm_perday$Variant_Name_NewNoPlus == "Delta")
cutdata_omicron<-subset(data_for_hm_perday, data_for_hm_perday$Variant_Name_NewNoPlus == "Omicron")

##trim data according to above dates when it goes below 50% consistently
cutdata_wt<-subset(cutdata_wt, cutdata_wt$day < "2021-03-09")
cutdata_alpha<-subset(cutdata_alpha, cutdata_alpha$day < "2021-06-20")
cutdata_delta<-subset(cutdata_delta, cutdata_delta$day < "2021-12-18")
cutdata_omicron<-subset(cutdata_omicron, cutdata_omicron$day < "2022-01-01")

##write out to file and edit to impute manually and to verify data is sensible 
##imputation reviewed in excel and performed as: mean of previous and next day which had real data values
###we decided skip WT due to lower amount of samples/runs from that time making the data gappy, and then start alpha at Feb 6 again.
#write.xlsx(cutdata_wt, "cutdata_dailyproportion_wt_trimmed.xlsx")
write.xlsx(cutdata_alpha, "cutdata_dailyproportion_alpha_trimmed.xlsx")
write.xlsx(cutdata_delta, "cutdata_dailyproportion_delta_trimmed.xlsx")
write.xlsx(cutdata_omicron, "cutdata_dailyproportion_omicron_trimmed.xlsx")

########alpha 5 day imputed########
cutdata_alpha<-read.xlsx("cutdata_dailyproportion_alpha_trimmed_imputed.xlsx")
cutdata_alpha$day<-as.Date(cutdata_alpha$Date.of.collection.Prcsd, origin = "1899-12-30")
#start alpha at Feb 6 since that's when coverage improves
cutdata_alpha<-cutdata_alpha[cutdata_alpha$day > "2021-02-05",]

#need to recalculate days span
cutdata_alpha$daysspan<-cutdata_alpha[1,6] %--% cutdata_alpha[,6]
cutdata_alpha$daysnum<-as.duration(cutdata_alpha$daysspan)/ ddays(1)

#set up table
alpha_processed_5days <- data.frame(StartDate=character(), EndDate=character(), StartDayNum=numeric(), EndDayNum=numeric(), 
                                            numDays=numeric(), Variant=character(), lmSlope=numeric(), lmIntercept=numeric())
nrow(cutdata_alpha) #134 rows
##use a 5 day window 
#small for loop to run basic linear model on the rolling window of data and save to a table with the intrcept 
for(i in c(seq(1, 129, 1))){
  startdate<-i
  enddate<-i+4
  datesubset<-cutdata_alpha[startdate:enddate,]
  estimate_datesubset <- lm((DailyPercent) ~ daysnum, data=datesubset) #measure days since point 1
  cf <- coef(estimate_datesubset)
  #cf will now contain a vector with the Intercept the slope
  estIntercept <- cf["(Intercept)"]
  estSlope <- cf[2]
  head(datesubset)
  alpha_processed_5days[nrow(alpha_processed_5days)+1,]<-
    c(as.character(datesubset[1,6]), as.character(datesubset[5,6]), datesubset[1,8], datesubset[5,8], (max(datesubset$daysnum) - min(datesubset$daysnum)), "alpha", estSlope, estIntercept)
    
} 

#now replot cutdata_alpha_v2 with x axis being days 0 to .. 
alpha_processed_5days$y1<-as.numeric(alpha_processed_5days$StartDayNum) * as.numeric(alpha_processed_5days$lmSlope) + 
      as.numeric(alpha_processed_5days$lmIntercept)
alpha_processed_5days$y2<-as.numeric(alpha_processed_5days$EndDayNum) * as.numeric(alpha_processed_5days$lmSlope) + 
      as.numeric(alpha_processed_5days$lmIntercept)

alpha_processed_5days[,c(3,9,4,10)]
head(alpha_processed_5days[,c(3,9,4,10)])
alpha_processed_5days$StartDayNum<-as.numeric(alpha_processed_5days$StartDayNum)
alpha_processed_5days$EndDayNum<-as.numeric(alpha_processed_5days$EndDayNum)

##save file
#write.xlsx(alpha_processed_5days, "finalized_alpha_processed_imputed_5days.xlsx")
alpha_processed_5days<-read.xlsx("finalized_alpha_processed_imputed_5days.xlsx")

alpha1<-ggplot(cutdata_alpha, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Alpha Cases, as Proportion of Cases Per Day") +
  xlab("Collection Day") + ylab("Daily Proportion") + theme_bw() + #scale_x_date(date_labels="%b %d %Y", date_breaks = "1 week") +
   geom_segment(data = alpha_processed_5days, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      #axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
alpha1
pdf("finalized_alpha_5daywindow_dailycases_withsegments_proportions_v2.pdf", width = 6, height = 4)
alpha1
dev.off()

##zoomed vers for better view of earlier days
cutdata_alpha_forzoomplot<-cutdata_alpha[cutdata_alpha$daysnum >= 0 & cutdata_alpha$daysnum <= 50,]
alpha_processed_5days_forzoomplot<-alpha_processed_5days[alpha_processed_5days$StartDate %in% cutdata_alpha_forzoomplot$day ,]

alpha1zoomNOLINE<-ggplot(cutdata_alpha_forzoomplot, aes(x=day, y = DailyPercent)) + geom_col() +
  ggtitle("Alpha Cases, as Proportion of Cases Per Day") +  ylim(0,1) +
  xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),    #axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8))
alpha1zoomNOLINE

#version with line showing the models
alpha1zoom<-ggplot(cutdata_alpha_forzoomplot, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Alpha Cases, as Proportion of Cases Per Day") + ylim(0,1) + 
  xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() +
   geom_segment(data = alpha_processed_5days_forzoomplot, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.y = element_text(hjust = 1, size = 8))
alpha1zoom

pdf("finalized_alpha_zoomed_sliding_window_5daywindow_imputed.pdf", width = 6, height = 4)
alpha1zoomNOLINE
alpha1zoom
dev.off()

#mini vers that only goes to when the variant is over 50% for 5 days.
#so for alpha it is going until March 6-March 10 2021
cutdata_alpha_forzoomplot50<-cutdata_alpha[cutdata_alpha$day <= "2021-03-10",]
alpha_processed_5days_forzoomplot<-alpha_processed_5days[alpha_processed_5days$StartDate %in% cutdata_alpha_forzoomplot50$day ,]

alpha1zoomto50percent<-ggplot(cutdata_alpha_forzoomplot50, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Alpha Cases, as Proportion of Cases Per Day") + ylim(0,1) + #scale_x_continuous(limits=c(-1,50)) + 
  xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() + #scale_x_date(date_labels="%b %d %Y", date_breaks = "1 week") +
   geom_segment(data = alpha_processed_5days_forzoomplot, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 14))
alpha1zoomto50percent

head(alpha_processed_5days_forzoomplot)
alpha_processed_5days_forzoomplot_50percent<-subset(alpha_processed_5days_forzoomplot, alpha_processed_5days_forzoomplot$StartDayNum %in% cutdata_alpha_forzoomplot50$daysnum)
alpha_processed_5days_forzoomplot_50percent$fillcol<-"positive"
alpha_processed_5days_forzoomplot_50percent[alpha_processed_5days_forzoomplot_50percent$lmSlope<0,]$fillcol<-"negative"
slopescols<-c("positive" = "red", negative = "blue")

alphachunkslopesDatesXZoom_50percent<-ggplot(alpha_processed_5days_forzoomplot_50percent, aes(x=as.Date(StartDate), y=as.numeric(lmSlope) )) + geom_col(aes(fill = fillcol)) +
  ggtitle("Alpha Window Slopes, Imputed Data") + xlab("") + ylab("Slope of Window") + theme_bw() + geom_line(size=2, alpha = 0.5) +
  scale_y_continuous(limits = c(-0.15,0.15)) + scale_fill_manual(values = slopescols) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),   axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 14)) 

pdf("finalized_alpha_5daywindow_zoomed_in_slopes_Imputed_Until50Percent_WithLine.pdf", width = 6, height = 4)
alpha1zoomto50percent
alphachunkslopesDatesXZoom_50percent
dev.off()

###full version of slopes for alpha 
alpha_processed_5days$fillcol<-"positive"
alpha_processed_5days[alpha_processed_5days$lmSlope<0,]$fillcol<-"negative"
slopescols<-c("positive" = "red", negative = "blue")

alphachunkslopesDates<-ggplot(alpha_processed_5days, aes(x=as.Date(StartDate), y=as.numeric(lmSlope), fill = fillcol ) ) + geom_col() + 
  scale_x_date(date_labels="%b %d %Y", date_breaks = "1 week") + scale_y_continuous(limits = c(-0.15,0.15)) +
  ggtitle("Alpha Window Slopes, Imputed Data") + xlab("Start Day") + ylab("Slope of Window") + theme_bw() + scale_fill_manual(values = slopescols) +
  theme(axis.title = element_text(size=10), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8)) 
pdf("finalized_alpha_5daywindow_slopes_Imputed_v2.pdf", width = 6, height = 4)
alphachunkslopesDates
dev.off()


####delta 5 day window, imputed######
cutdata_delta<-read.xlsx("cutdata_dailyproportion_delta_trimmed_imputed.xlsx")
cutdata_delta$day<-as.Date(cutdata_delta$Date.of.collection.Prcsd, origin = "1899-12-30")

#need to recalculate days span
cutdata_delta$daysspan<-cutdata_delta[1,6] %--% cutdata_delta[,6]
cutdata_delta$daysnum<-as.duration(cutdata_delta$daysspan)/ ddays(1)

#set up table
delta_processed_5days <- data.frame(StartDate=character(), EndDate=character(), StartDayNum=numeric(), EndDayNum=numeric(), 
                                            numDays=numeric(), Variant=character(), lmSlope=numeric(), lmIntercept=numeric())
nrow(cutdata_delta) #229 rows
##use a 5 day window 
for(i in c(seq(1, 225, 1))){
  startdate<-i
  enddate<-i+4
  datesubset<-cutdata_delta[startdate:enddate,]
  estimate_datesubset <- lm((DailyPercent) ~ daysnum, data=datesubset) 
  cf <- coef(estimate_datesubset)
  estIntercept <- cf["(Intercept)"]
  estSlope <- cf[2]
  head(datesubset)
  delta_processed_5days[nrow(delta_processed_5days)+1,]<-
    c(as.character(datesubset[1,6]), as.character(datesubset[5,6]), datesubset[1,8], datesubset[5,8], (max(datesubset$daysnum) - min(datesubset$daysnum)), "delta", estSlope, estIntercept)
    
} 

#organize and replot cutdata_delta_v2 with x axis being days 0 to ...
delta_processed_5days$y1<-as.numeric(delta_processed_5days$StartDayNum) * as.numeric(delta_processed_5days$lmSlope) + 
      as.numeric(delta_processed_5days$lmIntercept)
delta_processed_5days$y2<-as.numeric(delta_processed_5days$EndDayNum) * as.numeric(delta_processed_5days$lmSlope) + 
      as.numeric(delta_processed_5days$lmIntercept)

delta_processed_5days[,c(3,9,4,10)]
head(delta_processed_5days[,c(3,9,4,10)])
delta_processed_5days$StartDayNum<-as.numeric(delta_processed_5days$StartDayNum)
delta_processed_5days$EndDayNum<-as.numeric(delta_processed_5days$EndDayNum)

##save file
#write.xlsx(delta_processed_5days, "finalized_delta_processed_imputed_5days.xlsx")
delta_processed_5days<-read.xlsx("finalized_delta_processed_imputed_5days.xlsx")

#graph showing results of linear models
delta1<-ggplot(cutdata_delta, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Delta Cases, as Proportion of Cases Per Day") +  xlab("Collection Day") + ylab("Daily Proportion") + theme_bw() +
   geom_segment(data = delta_processed_5days, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.text.y = element_text(hjust = 1, size = 8))
delta1
pdf("finalized_delta_5daywindow_dailycases_withsegments_proportions_v2.pdf", width = 7, height = 4)
delta1
dev.off()

##zoomed vers to get better view of earlier days 
cutdata_delta_forzoomplot<-cutdata_delta[cutdata_delta$daysnum >= 0 & cutdata_delta$daysnum <= 50,]
delta_processed_5days_forzoomplot<-delta_processed_5days[delta_processed_5days$StartDate %in% cutdata_delta_forzoomplot$day ,]

delta1zoomNOLINE<-ggplot(cutdata_delta_forzoomplot, aes(x=day, y = DailyPercent)) + geom_col() +
  ggtitle("Delta Cases, as Proportion of Cases Per Day") +  ylim(0,1) +  xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() + 
   theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(hjust = 1, size = 8))
delta1zoomNOLINE

delta1zoom<-ggplot(cutdata_delta_forzoomplot, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Delta Cases, as Proportion of Cases Per Day") + ylim(0,1) +  xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() +
   geom_segment(data = delta_processed_5days_forzoomplot, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),   axis.text.y = element_text(hjust = 1, size = 8))
delta1zoom

pdf("finalized_delta_zoomed_sliding_window_5daywindow_imputed.pdf", width = 6, height = 4)
delta1zoomNOLINE
delta1zoom
dev.off()


#mini vers. that only goes to when the variant is over 50% for 5 days. so for delta it is going until June 16 to June 20 2021
cutdata_delta_forzoomplot50<-cutdata_delta[cutdata_delta$day <= "2021-06-20",]
delta_processed_5days_forzoomplot<-delta_processed_5days[delta_processed_5days$StartDate %in% cutdata_delta_forzoomplot50$day ,]

delta1zoomto50percent<-ggplot(cutdata_delta_forzoomplot50, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Delta Cases, as Proportion of Cases Per Day") + ylim(0,1) + xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() +
   geom_segment(data = delta_processed_5days_forzoomplot, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 14))
delta1zoomto50percent

head(delta_processed_5days_forzoomplot)
delta_processed_5days_forzoomplot_50percent<-subset(delta_processed_5days_forzoomplot, delta_processed_5days_forzoomplot$StartDayNum %in% cutdata_delta_forzoomplot50$daysnum)
delta_processed_5days_forzoomplot_50percent$fillcol<-"positive"
delta_processed_5days_forzoomplot_50percent[delta_processed_5days_forzoomplot_50percent$lmSlope<0,]$fillcol<-"negative"
slopescols<-c("positive" = "red", negative = "blue")

deltachunkslopesDatesXZoom_50percent<-ggplot(delta_processed_5days_forzoomplot_50percent, aes(x=as.Date(StartDate), y=as.numeric(lmSlope)) ) + geom_col(aes( fill = fillcol )) +
  ggtitle("Delta Window Slopes, Imputed Data") + xlab("") + ylab("Slope of Window") + theme_bw() +  scale_fill_manual(values = slopescols) +
   geom_line(size=2, alpha = 0.5) +  scale_y_continuous(limits = c(-0.15,0.15)) + theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 14)) 

pdf("finalized_delta_5daywindow_zoomed_in_slopes_Imputed_Until50Percent_withLine.pdf", width = 6, height = 4)
delta1zoomto50percent
deltachunkslopesDatesXZoom_50percent
dev.off()


###full version of slopes
delta_processed_5days$fillcol<-"positive"
delta_processed_5days[delta_processed_5days$lmSlope<0,]$fillcol<-"negative"
slopescols<-c("positive" = "red", negative = "blue")

deltachunkslopesDates<-ggplot(delta_processed_5days, aes(x=as.Date(StartDate), y=as.numeric(lmSlope), fill = fillcol ) ) + geom_col() + scale_fill_manual(values = slopescols) +
  ggtitle("Delta Window Slopes, Imputed Data") + xlab("Start Day") + ylab("Slope of Window") + theme_bw() + scale_x_date(date_labels="%b %d %Y", date_breaks = "1 week") +
  scale_y_continuous(limits = c(-0.15,0.15)) + theme(axis.title = element_text(size=10), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8)) 
pdf("finalized_delta_5daywindow_slopes_Imputed_v2.pdf", width = 7, height = 4)
deltachunkslopesDates
dev.off()


#####omicron 5 day window, imputed#####
cutdata_omicron<-read.xlsx("cutdata_dailyproportion_omicron_trimmed_imputed.xlsx")
cutdata_omicron$day<-as.Date(cutdata_omicron$Date.of.collection.Prcsd, origin = "1899-12-30")

#need to recalculate days span
cutdata_omicron$daysspan<-cutdata_omicron[1,6] %--% cutdata_omicron[,6]
cutdata_omicron$daysnum<-as.duration(cutdata_omicron$daysspan)/ ddays(1)

#set up table
omicron_processed_5days <- data.frame(StartDate=character(), EndDate=character(), StartDayNum=numeric(), EndDayNum=numeric(), 
                                            numDays=numeric(), Variant=character(), lmSlope=numeric(), lmIntercept=numeric())
nrow(cutdata_omicron) #33 rows
##5 day window for analysis. shorter here because we end the omicron analysis around Dec 31 because that is when the government ended testing for the public 
for(i in c(seq(1, 29, 1))){
  startdate<-i
  enddate<-i+4
  datesubset<-cutdata_omicron[startdate:enddate,]
  estimate_datesubset <- lm((DailyPercent) ~ daysnum, data=datesubset) 
  cf <- coef(estimate_datesubset)
  estIntercept <- cf["(Intercept)"]
  estSlope <- cf[2]
  head(datesubset)
  omicron_processed_5days[nrow(omicron_processed_5days)+1,]<-
    c(as.character(datesubset[1,6]), as.character(datesubset[5,6]), datesubset[1,8], datesubset[5,8], (max(datesubset$daysnum) - min(datesubset$daysnum)), "omicron", estSlope, estIntercept)
} 

#now replot cutdata_omicron_v2 with x axis being days 0 to ...
omicron_processed_5days$y1<-as.numeric(omicron_processed_5days$StartDayNum) * as.numeric(omicron_processed_5days$lmSlope) + 
      as.numeric(omicron_processed_5days$lmIntercept)
omicron_processed_5days$y2<-as.numeric(omicron_processed_5days$EndDayNum) * as.numeric(omicron_processed_5days$lmSlope) + 
      as.numeric(omicron_processed_5days$lmIntercept)

omicron_processed_5days[,c(3,9,4,10)]
head(omicron_processed_5days[,c(3,9,4,10)])
omicron_processed_5days$StartDayNum<-as.numeric(omicron_processed_5days$StartDayNum)
omicron_processed_5days$EndDayNum<-as.numeric(omicron_processed_5days$EndDayNum)

##save file
write.xlsx(omicron_processed_5days, "finalized_omicron_processed_imputed_5days.xlsx")

omicron1<-ggplot(cutdata_omicron, aes(x=daysnum, y = DailyPercent)) + geom_col() +  ggtitle("Omicron") +
  xlab("Collection Day") + ylab("Daily Proportion") + theme_bw() + 
   geom_segment(data = omicron_processed_5days, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.y = element_text(hjust = 1, size = 8))
omicron1
pdf("finalized_omicron_5daywindow_dailycases_withsegments_proportions_v2.pdf", width = 2, height = 4)
omicron1
dev.off()

##zoomed vers for better view of first few days - actually the same for omicron as it is a short time period before testing was cut 
cutdata_omicron_forzoomplot<-cutdata_omicron[cutdata_omicron$daysnum >= 0 & cutdata_omicron$daysnum <= 50,]
omicron_processed_5days_forzoomplot<-omicron_processed_5days[omicron_processed_5days$StartDate %in% cutdata_omicron_forzoomplot$day ,]

omicron1zoomNOLINE<-ggplot(cutdata_omicron_forzoomplot, aes(x=day, y = DailyPercent)) + geom_col() +
  ggtitle("Omicron Cases, as Proportion of Cases Per Day") +  ylim(0,1) +  xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(hjust = 1, size = 8))
omicron1zoomNOLINE

omicron1zoom<-ggplot(cutdata_omicron_forzoomplot, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Omicron Cases, as Proportion of Cases Per Day") + ylim(0,1) + xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() + 
   geom_segment(data = omicron_processed_5days_forzoomplot, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(hjust = 1, size = 8))
omicron1zoom

pdf("finalized_omicron_zoomed_sliding_window_5daywindow_imputed.pdf", width = 6, height = 4)
omicron1zoomNOLINE
omicron1zoom
dev.off()


#new mini vers that only goes to when the variant is over 50% for 5 days; so for omicron it is going until Dec 14-18 2021
cutdata_omicron_forzoomplot50<-cutdata_omicron[cutdata_omicron$day <= "2021-12-18",]
omicron_processed_5days_forzoomplot<-omicron_processed_5days[omicron_processed_5days$StartDate %in% cutdata_omicron_forzoomplot50$day ,]

omicron1zoomto50percent<-ggplot(cutdata_omicron_forzoomplot50, aes(x=daysnum, y = DailyPercent)) + geom_col() +
  ggtitle("Omicron Cases, as Proportion of Cases Per Day") + ylim(0,1) + xlab("Collection Day") + ylab("Daily Proportion of Cases") + theme_bw() +
   geom_segment(data = omicron_processed_5days_forzoomplot, aes(x = StartDayNum, y = y1, xend = EndDayNum, yend = y2), colour="red") +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 14))
omicron1zoomto50percent

head(omicron_processed_5days_forzoomplot)
omicron_processed_5days_forzoomplot_50percent<-subset(omicron_processed_5days_forzoomplot, omicron_processed_5days_forzoomplot$StartDayNum %in% cutdata_omicron_forzoomplot50$daysnum)
omicron_processed_5days_forzoomplot_50percent$fillcol<-"positive"
#omicron_processed_5days_forzoomplot_50percent[omicron_processed_5days_forzoomplot_50percent$lmSlope<0,]$fillcol<-"negative" #no negative in omicron 
slopescols<-c("positive" = "red", negative = "blue")

omicronchunkslopesDatesXZoom_50percent<-ggplot(omicron_processed_5days_forzoomplot_50percent, aes(x=as.Date(StartDate), y=as.numeric(lmSlope)) ) + geom_col(aes( fill = fillcol )) +
  ggtitle("Omicron Window Slopes, Imputed Data") + xlab("") + ylab("Slope of Window") + theme_bw() +  scale_fill_manual(values = slopescols) +
  scale_y_continuous(limits = c(-0.15,0.15)) + geom_line(size=2, alpha = 0.5) +
  theme(axis.title = element_text(size=12), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),   axis.text.x = element_blank(),axis.text.y = element_text(hjust = 1, size = 14)) 

pdf("finalized_omicron_5daywindow_zoomed_in_slopes_Imputed_Until50Percent_withLine.pdf", width = 6, height = 4)
omicron1zoomto50percent
omicronchunkslopesDatesXZoom_50percent
dev.off()

###full version of slopes
omicron_processed_5days$fillcol<-"positive"
omicron_processed_5days[omicron_processed_5days$lmSlope<0,]$fillcol<-"negative"
slopescols<-c("positive" = "red", negative = "blue")


omicronchunkslopesDates<-ggplot(omicron_processed_5days, aes(x=as.Date(StartDate), y=as.numeric(lmSlope), fill = fillcol ) ) + geom_col() + scale_fill_manual(values = slopescols) +
  ggtitle("Omicron") + xlab("Start Day") + ylab("Slope of Window") + theme_bw() + scale_x_date(date_labels="%b %d %Y", date_breaks = "1 week") + scale_y_continuous(limits = c(-0.15,0.15)) +
  theme(axis.title = element_text(size=10), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 8)) 
pdf("finalized_omicron_5daywindow_slopes_Imputed_v2.pdf", width = 2, height = 4)
omicronchunkslopesDates
dev.off()

#####end of alpha/delta/omicron individual analysis#####

####next, convert the slopes into a heatmap#######
###merge the 3 tables
all3_5days<-rbind(alpha_processed_5days, delta_processed_5days, omicron_processed_5days)
write.xlsx(all3_5days, file="alpha_delta_omicron_processed_5daywindows_imputed.xlsx")
all3_5days<-read.xlsx("alpha_delta_omicron_processed_5daywindows_imputed.xlsx")
all3_5days$Variant_F<-factor(all3_5days$Variant, levels=c("omicron", "delta", "alpha"))

all3_heatmap<-ggplot(all3_5days, aes(x=as.Date(StartDate), y=Variant_F, fill = as.numeric(lmSlope) )) + geom_tile() +
  ggtitle("Alpha, Delta, Omicron 5-Day Window Slopes, incl. Imputed Data") + xlab("Start Day of Window") + ylab("Variant") + theme_bw() + 
  scale_fill_gradientn(name = "Slope", colours = c("blue","gray","red")) + scale_x_date(date_labels="%b %d %Y", date_breaks = "2 weeks") +
  ##add two lines - one marking start of daily sparseq  at June 1 2021 and one marking end of mass testing at Dec 31 2021
  geom_vline(xintercept = as.Date("2021-06-01"), linetype = "dashed" ) +
  geom_vline(xintercept = as.Date("2021-12-31"), linetype = "dashed", colour = "red") +
  theme(axis.title = element_text(size=12), legend.position = c("bottom"), panel.background = element_rect(fill="white"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(angle = 45, hjust=1),
      plot.title = element_text(hjust=0.5),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  axis.text.y = element_text(hjust = 1, size = 12)) 

pdf("finalized_all3_5daywindow_slopes_Imputed_heatmap_v4_revisionsVers.pdf", width = 12, height = 4)
all3_heatmap
dev.off()



####end of heatmap#####
###end of finalized vers of figure 3 analysis####
