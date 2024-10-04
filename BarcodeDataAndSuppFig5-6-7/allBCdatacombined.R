

library("ggplot2")
library("dplyr")
library("openxlsx")
library("stringr")


###begin with stat tables from individual processing of barcode results
#then we can look at the contamination and spread stats per week 

run27stats<-read.xlsx("PythonBC_stats_Run27.xlsx", sheet = "Correct BCs in Correct Well %")
run27spread<-read.xlsx("PythonBC_stats_Run27.xlsx", sheet = "Correct BCs in Wrong Well %")
#edit triple sample
run27stats[run27stats$sampleID == "W60603827",]$sampleID<-paste("W60603827", run27stats[run27stats$sampleID == "W60603827",]$BC_Well, sep = "_")

run31stats<-read.xlsx("PythonBC_stats_Run31.xlsx", sheet = "Correct BCs in Correct Well %")
run31spread<-read.xlsx("PythonBC_stats_Run31.xlsx", sheet = "Correct BCs in Wrong Well %")

run38stats<-read.xlsx("PythonBC_stats_Run38.xlsx", sheet = "Correct BCs in Correct Well %")
run38spread<-read.xlsx("PythonBC_stats_Run38.xlsx", sheet = "Correct BCs in Wrong Well %")

run46stats<-read.xlsx("PythonBC_stats_Run46.xlsx", sheet = "Correct BCs in Correct Well %")
run46spread<-read.xlsx("PythonBC_stats_Run46.xlsx", sheet = "Correct BCs in Wrong Well %")

run139stats<-read.xlsx("PythonBC_stats_139.xlsx", sheet = "Correct BCs in Correct Well %")
run139spread<-read.xlsx("PythonBC_stats_139.xlsx", sheet = "Correct BCs in Wrong Well %")

run141p3stats<-read.xlsx("PythonBC_stats_141P3.xlsx", sheet = "Correct BCs in Correct Well %")
run141p3spread<-read.xlsx("PythonBC_stats_141P3.xlsx", sheet = "Correct BCs in Wrong Well %")

run141p6stats<-read.xlsx("PythonBC_stats_141P6.xlsx", sheet = "Correct BCs in Correct Well %")
run141p6spread<-read.xlsx("PythonBC_stats_141P6.xlsx", sheet = "Correct BCs in Wrong Well %")

#for 139 and 141 need to drop well after _
run141p3stats$sampleID<-gsub("_.*","",run141p3stats$sampleID)
run139stats$sampleID<-gsub("_.*","",run139stats$sampleID)
run141p6stats$sampleID<-gsub("_.*","",run141p6stats$sampleID)

run27stats<-run27stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run31stats<-run31stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run38stats<-run38stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run46stats<-run46stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run139stats<-run139stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run141p3stats<-run141p3stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run141p6stats<-run141p6stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
 
#drop controls
allstats<-rbind(run27stats, run31stats, run38stats, run46stats, run139stats, run141p3stats, run141p6stats)
allstats<-allstats[grep("^H", allstats$sampleID, invert = T),]
allstats<-allstats[grep("^N", allstats$sampleID, invert = T),]
allstats<-allstats[grep("^X00000000", allstats$sampleID, invert = T),]

##load manual annotations
###note that here anything that isn't the main 4 (wt/alpha/delta/omi) is called Other. so includes ind., mix, etc
annos<-read.xlsx("BC_test_annotations_for_finalfigure.xlsx")
allstats<-merge(allstats, annos, by.x = "sampleID", by.y = "sample", all = T)
#drop dupe, excluded, negative
allstats<-subset(allstats, !(allstats$BC_call %in% c("dupe", "excluded", "negative")))

allstats$run_c<-as.character(allstats$run)
allstats$run_F<-factor(allstats$run_c, levels =c("27", "31", "38", "46", "139", "141"))

allstats[is.na(allstats$percent_of_counts_within_sample),]
#remove the few samples which have NA for python scores
# wells that had 0 of their own BCs so they don't even have a score at all 
#this would be the ones who are blank on the heatmaps
#so I should set them to show as 0, not to be NA and remove them
allstats[is.na(allstats$percent_of_counts_within_sample),]$percent_of_counts_within_sample<-0

###add inverse which is the amount of intrawell contamination so other BCs that came into this well
allstats$intrawell<-(100-allstats$percent_of_counts_within_sample)
#save
write.xlsx(allstats, "allstatsforCorrectBCinCorrectWell_allBCruns_v1.xlsx")


###repeat this  for spread####
run139spread$Experiment<-139
run141p3spread$Experiment<-141
run141p6spread$Experiment<-141
run27spread$Experiment<-27
run31spread$Experiment<-31
run38spread$Experiment<-38
run46spread$Experiment<-46

#merge with other data 
#merge by BC_well and run 
run27spread<-merge(run27spread, run27stats, by ="BC_Well", all=T)
run31spread<-merge(run31spread, run31stats, by ="BC_Well", all=T)
run38spread<-merge(run38spread, run38stats, by ="BC_Well", all=T)
run46spread<-merge(run46spread, run46stats, by ="BC_Well", all=T)
run139spread<-merge(run139spread, run139stats, by ="BC_Well", all=T)
run141p3spread<-merge(run141p3spread, run141p3stats, by ="BC_Well", all=T)
run141p6spread<-merge(run141p6spread, run141p6stats, by ="BC_Well", all=T)
allspread<-rbind(run27spread, run31spread, run38spread, run46spread, run139spread, run141p3spread, run141p6spread)
allspread$percent_of_counts_within_sample<-NULL 

annos<-read.xlsx("BC_test_annotations_for_finalfigure.xlsx")

allspread<-allspread[grep("^H", allspread$sampleID, invert = T),]
allspread<-allspread[grep("^N", allspread$sampleID, invert = T),]
allspread<-allspread[grep("^X00000000", allspread$sampleID, invert = T),]

allspread<-merge(allspread, annos, by.x = "sampleID", by.y = "sample", all = T)
#drop dupe, excluded, negative
allspread<-subset(allspread, !(allspread$BC_call %in% c("dupe", "excluded", "negative")))
table(allspread$BC_call)


allspread$run_c<-as.character(allspread$Experiment)
allspread$run_F<-factor(allspread$run_c, levels =c("27", "31", "38", "46", "139", "141"))
allspread<-subset(allspread, !(is.na(allspread$percentTotal)))
allspread[is.na(allspread$BC_call),]

#     sampleID BC_Well percentTotal Experiment Variant_Name BC_call run run_c run_F
#2002     <NA>     11A          100         38         <NA>     <NA>  NA    38    38 >> W72902322 bad sample, alpha
#2003     <NA>     23K          100         46         <NA>     <NA>  NA    46    46 >> control
#2004     <NA>      9J          100        141         <NA>     <NA>  NA   141   141 >> X61609238 bad sample, omicron
#2005     <NA>     23E          100         38         <NA>     <NA>  NA    38    38 >> control
#2006     <NA>     23G          100         38         <NA>     <NA>  NA    38    38 >> control
#2007     <NA>     22G          100         38         <NA>     <NA>  NA    38    38 >> W73008922 bad sample, alpha 

#so fix the fixable wells and then remove the controls 
allspread[allspread$Experiment == 38 & allspread$BC_Well == "11A",]$BC_call<-"Alpha"
allspread[allspread$Experiment == 141 & allspread$BC_Well == "9J",]$BC_call<-"Omicron" #oddly only picks up actual 9J not other 141 plate 9J
allspread[allspread$Experiment == 38 & allspread$BC_Well == "22G",]$BC_call<-"Alpha"

#drop controls
allspreadnocontrols<-subset(allspread, !(allspread$BC_Well %in% c("23K", "23E", "23G")))

#so this is amount of a well;s BCs that travelled around. 
#save
write.xlsx(allspread, "BC1_allspreadforCorrectBCinIncorrectWell_allBCruns_v1.xlsx")

###make merge of stat and spread tables for sup table 1##
head(allstats)
head(allspread)

combodata<-merge(allstats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample", "intrawell", "Variant_Name",  "BC_call", "run")], 
                 allspread[,c("sampleID", "percentTotal", "Experiment" )], by="sampleID")
head(combodata)
#sampleID	BC_Well	percent_of_counts_within_sample	intrawell	Variant_Name	BC_call	run	percentTotal	Experiment
colnames(combodata)<-c("sampleID",	"BC_Well",	"percent_of_counts_within_sample",	"intrawell_contamination",	"Variant_Name",	"BC_call",	"run",	"interwell_contamination",	"Experiment")
write.xlsx(combodata, file="combined_spread_stats_nocontrols.xlsx")
##this is basically supplemental data 1 just with more columns # 


#####boxplots with scatters for supp. fig 5#####

##load annos
###note that here anything that isn't the main 4 (wt/alpha/delta/omi) is called Other. so includes ind., mix, etc
annos<-read.xlsx("BC_test_annotations_for_finalfigure.xlsx")

#load table
allstats<-read.xlsx("allstatsforCorrectBCinCorrectWell_allBCruns_v1.xlsx")
allstats$run_F<-factor(allstats$run_F, levels = c("27", "31", "38", "46", "139", "141"))

###version with change beta/gamma, Indeterminate, Mu, to other
allstatsv2[allstatsv2$BC_call %in% c("Beta/Gamma", "Indeterminate", "Mu"),]$BC_call<-"Other"
allstatsv2[allstatsv2$percent_of_counts_within_sample == 0,]


pdf("BCdata_correctlypaired_Boxplots_byplate.pdf", width = 6, height = 5)
python_allwells_correctPair_violin2<-ggplot(allstatsv2, aes(x = run_F, y=percent_of_counts_within_sample)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.2, aes(colour = BC_call)) + 
  ggtitle("BC Data:\nPaired R1+R2 Results\n% of Correctly Paired BCs in Correct Well (per well)") + 
  xlab("") + ylab("Percent") + theme_bw() +  ylim(0,100) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10, hjust = 0.5),  axis.text.y = element_text(hjust = 1, size = 10))
python_allwells_correctPair_violin2
dev.off()

  
####repeat this graph for spread####
allspread<-read.xlsx("BC1_allspreadforCorrectBCinIncorrectWell_allBCruns_v1.xlsx")
allspread$run_F<-factor(allspread$run_F, levels = c("27", "31", "38", "46", "139", "141"))

allspreadv2<-allspread
allspreadv2[allspreadv2$BC_call %in% c("Beta/Gamma", "Indeterminate", "Mu"),]$BC_call<-"Other"
allspreadv2[is.na(allspreadv2$BC_call),]
allspreadv2<-subset(allspreadv2,!(allspreadv2$BC_Well %in% c("23K", "23E", "23G")))

pdf("BCdata_correctlypairedIncorrectWell_Boxplots_byplate_MoreOthers_v3.pdf", width = 6, height = 5)
spread_violin2<-ggplot(allspreadv2, aes(x = run_F, y=percentTotal)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.2, aes(colour = BC_call)) + 
  ggtitle("BC Data:\nPaired R1+R2 Results\n% of Correctly Paired BCs in Incorrect Well\n(per BC pair)") +
  xlab("") + ylab("Percent") + theme_bw() + 
  theme(axis.title = element_text(size=12), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10, hjust = 0.5),   axis.text.y = element_text(hjust = 1, size = 10))
spread_violin2
dev.off()


######log bowtie plot#####
#if we plot log(bowtiecounts+1) for each well and rank from highest to lowest on the x-axis it would be great 
###visually to show that the outliers of the violin plot are all clustered to the super low end

count27<-read.xlsx("run27_CountTable.xlsx"); count27$filename<-NULL
count31<-read.xlsx("run31_CountTable.xlsx"); count31$filename<-NULL
count38<-read.xlsx("run38_CountTable.xlsx"); count38$filename<-NULL
count46<-read.xlsx("run46_CountTable.xlsx"); count46$filename<-NULL
count139<-read.xlsx("runcl139_CountTable.xlsx"); count139$filename<-NULL
count141p3<-read.xlsx("runcl141P3_CountTable.xlsx"); count141p3$filename<-NULL
count141p6<-read.xlsx("runcl141P6_CountTable.xlsx"); count141p6$filename<-NULL

count27$run<-27
count31$run<-31
count38$run<-38
count46$run<-46
count139$run<-139
count141p3$run<-141
count141p6$run<-141

allcount<-rbind(count27, count31, count38, count46, count139, count141p3, count141p6)
allcount$sample<-gsub("-V1-2", "", allcount$sample)

#drop controls and negs
allcount<-allcount[grep("^H", allcount$sample, invert = T),]
allcount<-allcount[grep("^N", allcount$sample, invert = T),]
allcount<-allcount[grep("^X00000000", allcount$sample, invert = T),]
head(allcount)
tail(allcount, 20)

allcountplot1<-allcount
allcountplot1$well
allcountplot1$log10bowtie<-log10(allcountplot1$total.raw.reads+1)
allcountplot1 <- allcountplot1[order(allcountplot1$log10bowtie),]
allcountplot1$ordered_sample<-factor(allcountplot1$sample, levels=c(allcountplot1$sample))
allcountplot1$run<-as.character(allcountplot1$run)

write.xlsx(allcountplot1, "allBC1_log10bowtiecountlist.xlsx")

#merge with annos again
annos<-read.xlsx("../BC_test_annotations_for_finalfigure.xlsx"); annos$run<-NULL
allcount<-merge(allcount, annos, by.x = "sample", by.y = "sample", all = T)
#drop dupe, excluded, negative
allcount<-subset(allcount, !(allcount$BC_call %in% c("dupe", "excluded", "negative")))
table(allcount$BC_call)
#        Alpha    Beta/Gamma         Delta Indeterminate            Mu       Omicron         Other            WT 
#          585            50           205           148             4           456           240           292 

allcount$run_c<-as.character(allcount$run)
table(allcount$run)
allcount$run_F<-factor(allcount$run_c, levels =c("27", "31", "38", "46", "139", "141"))


allcount$log10bowtie<-log10(allcount$total.raw.reads+1)
allcount <- allcount[order(allcount$log10bowtie),]

#write list 
write.xlsx(allcount, "all_BC1_log10bowtiecountlist_withannotation.xlsx")
#list of red ones
redlistind_alldata<-subset(allcount, allcount$log10bowtie < 4.5)
write.xlsx(redlistind_alldata, file = "listofreddotsonlogcountplot_alldata.xlsx")


######supp fig 6######
###make all the log10 count plots 

allcountplot1<-read.xlsx("all_BC1_log10bowtiecountlist_withannotation.xlsx")
allcountplot1$well
allcountplot1$log10bowtie<-log10(allcountplot1$total.raw.reads+1)
allcountplot1 <- allcountplot1[order(allcountplot1$log10bowtie),]
allcountplot1$ordered_sample<-factor(allcountplot1$sample, levels=c(allcountplot1$sample))
allcountplot1$run<-as.character(allcountplot1$run)
allcountplot1$run_f<-factor(allcountplot1$run, levels=c("27", "31", "38", "46", "139", "141" ))
head(allcountplot1)

plotrun_orderedwells<-ggplot(allcountplot1, aes(x=ordered_sample, y=log10bowtie)) + facet_wrap(.~run_f, ncol=1) +
  geom_point(alpha = 0.75, size = 0.5) + ggtitle("BC Data, Bowtie Counts\nLog10(Total Raw Reads + 1)") + 
  xlab("Sample") + ylab("Log10 (Total Raw Reads + 1)") + theme_bw() + #scale_x_log10(breaks = c(1,10,100,1000,10000,100000, 1000000)) + 
  labs(colour="Run ID", shape = "Data Type") +
  geom_point(data = allcountplot1[allcountplot1$log10bowtie < 4.5,], color = "red", size = 0.5) +
  theme(axis.title = element_text(size=16), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(hjust = 1, size = 10))
plotrun_orderedwells


pdf("allBC_log10bowtiecounts_ordered.pdf", width = 10, height = 8)
plotrun_orderedwells
dev.off()

###end supp fig 6####



####supp fig 7: outlined heatmaps######
setwd("/Volumes/PromisePegasus/SPAR_SEQ/revisions_Sep2024/bc_processing_and_supfig567/bc1")

###begin with stat tables 
#then we can look at the contamination and spread stats per wekk 
run27stats<-read.xlsx("PythonBC_stats_Run27.xlsx", sheet = "Correct BCs in Correct Well %")
run27stats[run27stats$sampleID == "W60603827",]$sampleID<-paste("W60603827", run27stats[run27stats$sampleID == "W60603827",]$BC_Well, sep = "_")
run31stats<-read.xlsx("PythonBC_stats_Run31.xlsx", sheet = "Correct BCs in Correct Well %")
run38stats<-read.xlsx("PythonBC_stats_Run38.xlsx", sheet = "Correct BCs in Correct Well %")
run46stats<-read.xlsx("PythonBC_stats_Run46.xlsx", sheet = "Correct BCs in Correct Well %")
run139stats<-read.xlsx("PythonBC_stats_139.xlsx", sheet = "Correct BCs in Correct Well %")
run141p3stats<-read.xlsx("PythonBC_stats_141P3.xlsx", sheet = "Correct BCs in Correct Well %")
run141p6stats<-read.xlsx("PythonBC_stats_141P6.xlsx", sheet = "Correct BCs in Correct Well %")
#for 139 and 141 need to drop well after _
run141p3stats$sampleID<-gsub("_.*","",run141p3stats$sampleID)
run139stats$sampleID<-gsub("_.*","",run139stats$sampleID)
run141p6stats$sampleID<-gsub("_.*","",run141p6stats$sampleID)

run27stats<-run27stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run31stats<-run31stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run38stats<-run38stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run46stats<-run46stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run139stats<-run139stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run141p3stats<-run141p3stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]
run141p6stats<-run141p6stats[,c("sampleID", "BC_Well", "percent_of_counts_within_sample")]

#load run info since the other info is for samples only, not controls
run27stats$runID<-"27"
run31stats$runID<-"31"
run38stats$runID<-"38"
run46stats$runID<-"46"
run139stats$runID<-"139"
run141p3stats$runID<-"141p3"
run141p6stats$runID<-"141p6"

allstats<-rbind(run27stats, run31stats, run38stats, run46stats, run139stats, run141p3stats, run141p6stats)


##load annos
###note that here anything that isn't the main 4 (wt/alpha/delta/omi) is called Other. so includes ind., mix, etc
annos<-read.xlsx("../BC_test_annotations_for_finalfigure.xlsx")
allstats<-merge(allstats, annos, by.x = "sampleID", by.y = "sample", all = T)
allstats$run_c<-as.character(allstats$run)

allstats[is.na(allstats$percent_of_counts_within_sample),]
#remove the few samples which have NA for python scores
#they are wells that had 0 of their own BCs so they don't even have a score at all 
#this would be the ones who are blank on the heatmaps
#so I should set them to show as 0, not to be NA and remove them
allstats[is.na(allstats$percent_of_counts_within_sample),]$percent_of_counts_within_sample<-0
#finally fix runID as needed
allstats[is.na(allstats$runID),]$runID<-allstats[is.na(allstats$runID),]$run
#X61609238 is 141p3
allstats[allstats$sampleID == "X61609238",]$runID<-"141p3" 
allstats$run<-NULL
allstats$run_c<-NULL

#locate missing well ID for the samples that are NA/0%
#W72902322  11A        Alpha    Alpha  38    38
#W73008922  22G        Alpha    Alpha  38    38
#X61407567  21N       Delta    Delta 139   139
#X61508862  2J   Omicron    Delta 139   139
#X61508906  4L    Omicron    Delta 139   139
#X61609238  9J   141p3

allstats[allstats$sampleID == "W72902322", ]$BC_Well<-"11A"
allstats[allstats$sampleID == "W73008922", ]$BC_Well<-"22G"
allstats[allstats$sampleID == "X61407567", ]$BC_Well<-"21N"
allstats[allstats$sampleID == "X61508862", ]$BC_Well<-"2J"
allstats[allstats$sampleID == "X61508906", ]$BC_Well<-"4L"
allstats[allstats$sampleID == "X61609238", ]$BC_Well<-"9J"


#save
write.xlsx(allstats, "allstats_supFig7.xlsx")
#load low count samples list
redlistind_alldata<-read.xlsx("listofreddotsonlogcountplot_alldata.xlsx")

##heatmap with violin plot values
#colour border of wells for the ones that are red on the violin plot

#use pch=21, then colour is border and fill is middle colour
colourscale<-c("Fail" = "red", "Pass" = "black")

head(allstats)
allstats$percent<-allstats$percent_of_counts_within_sample

####This is the version in Sup Fig 7####
pdf("AllBC1Data_SuppFig7Heatmaps.pdf", width = 7, height = 5)
for (i in c("27", "31", "38", "46", "139", "141p3", "141p6")) {
  print(i)
  currentexp<-subset(allstats, allstats$runID == i)
  
  currentexp$row<-gsub("[0-9]", "", currentexp$BC_Well)
  currentexp$col<-gsub("[A-Z]", "", currentexp$BC_Well)
  currentexp$ordered_row<-factor(currentexp$row, levels=LETTERS[seq( from = 1, to = 16 )])
  currentexp$ordered_col<-factor(currentexp$col, levels=c(1:24))
  currentexp$BCcountscutoff<-"Pass"
  #we don't have controls on the red list 
  if (any(currentexp$sampleID %in% redlistind_alldata$sample )) {
    currentexp[currentexp$sampleID %in% redlistind_alldata$sample, ]$BCcountscutoff<-"Fail"}
  
  python_x_hm1<-ggplot(currentexp, aes(x=ordered_col, y=ordered_row, fill = percent, colour = BCcountscutoff)) +
    scale_fill_gradient(low = "blue", high = "white", na.value = "grey", breaks =c(0,25,50,75,100), limits = c(0,100)) + scale_y_discrete(limits=rev, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    geom_point(size = 4, pch=21, stroke=2) + ggtitle(paste("BC Run", i, ": % of Correctly Paired BCs in Correct Well, Per Well")) + 
    xlab("Column") + ylab("Row") + theme_bw() + 
    guides(colourbar = guide_legend(title.position = "top")) + scale_colour_manual(values = colourscale) +
    theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(hjust = 1, size = 10))
  print(python_x_hm1)
}
dev.off()

###end of Sup Fig 7#####


