######Supp Fig 8, large gray heatmap showing every variant########
#ordered over time
alldatareload<-read.xlsx("filtered_data_for_paper_v25.xlsx")
# 
alldatareload$Date.of.collection<-as.Date(alldatareload$Date.of.collection, origin = "1899-12-30") #
alldatareload<-alldatareload[,c("sample", "Variant.Details", "Date", "Date.of.collection",  "Variant_Name_New")]

#now need to use subvar list from above
#this is from Fig 5 script. 
subvarlist_forplot<-read.xlsx("organized_subvariants_v25.xlsx")

#drop main variants.
mainvars<-c("B.1.1.7", "B.1.1.529", "B.1.1.529 (BA.2/BA.3)", "B.1.1.529 (BA.4/BA.5)", "WT", "XBB", 
            "XBB.1.5", "B.1.1.525", "B.1.621", "B.1.617.2", "P1/B.1.351",
            "B.1.1.529 (BA.2.75)", "B.1.1.529 (BA.2.75.2)")

subvarlist_forplot<-subset(subvarlist_forplot, !(subvarlist_forplot$Variant.Details %in% mainvars))

submainlist<-subset(alldatareload, alldatareload$Variant.Details %in% subvarlist_forplot$Variant.Details)

#remove vars from variant details
subvarlist_forplot<-subvarlist_forplot
subvarlist_forplot$Variant.Details2<-subvarlist_forplot$Variant.Details
subvarlist_forplot$Variant.Details2<-gsub("B.1.1.529 \\(BA.4/BA.5\\)\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("B.1.1.529 \\(BA.2.75\\)\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("B.1.1.7\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("P1/B.1.351\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("B.1.617.2\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("B.1.1.529\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("B.1.1.529 \\(BA.2.75.2\\)\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("B.1.1.529 \\(BA.2/BA.3\\)\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("XBB.1.5\\+","",subvarlist_forplot$Variant.Details2)
subvarlist_forplot$Variant.Details2<-gsub("XBB\\+","",subvarlist_forplot$Variant.Details2)
head(subvarlist_forplot, 80)
subvarlist_forplot[subvarlist_forplot$Variant.Details2 == "",]
subvarlist_forplot<-subvarlist_forplot[subvarlist_forplot$Variant.Details2!="",]

table(subvarlist_forplot$Variant_Name_New_Cut)
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Other",]
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Omicron BA.2.75.2",]
###edge case - PBS insertion under Variant_Name_New_Cut = Other as well as B.1.617.1
##so my ordering for variant main is actually:
# WT, Other, Alpha, Beta/Gamma, Delta, Omicron, Omicron BA.2.75, Omicron BA.2.75.2, Omicron BA.2/3,  Omicron BA.4/5,  XBB, XBB.1.5

#add value for order position so WT = 1, Other = 2, ...
subvarlist_forplot$ordervec<-1
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Other",]$ordervec<-2
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Alpha",]$ordervec<-3
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Beta/Gamma",]$ordervec<-4
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Delta",]$ordervec<-5
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Eta",]$ordervec<-6
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Omicron",]$ordervec<-7
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Omicron BA.2.75",]$ordervec<-9
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Omicron BA.2/3",]$ordervec<-8
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "Omicron BA.4/5",]$ordervec<-10 #9
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "XBB",]$ordervec<-11 #10
subvarlist_forplot[subvarlist_forplot$Variant_Name_New_Cut == "XBB.1.5",]$ordervec<-12 #11


##next get order of pos numbers (first pos only in case of pairs)
subvarlist_forplot$firstpos<-gsub("[A-Za-z]", "", subvarlist_forplot$Variant.Details2)
#then take the first 3 chars
subvarlist_forplot$firstpos<-substr(subvarlist_forplot$firstpos, start = 1, stop = 3)
#order by #1 main var and then #2 position 
subvarlist_forplot<-subvarlist_forplot[order(subvarlist_forplot$ordervec, subvarlist_forplot$firstpos, decreasing = FALSE),]
##combine main var plus the variant info 
head(subvarlist_forplot)
subvarlist_forplot$MainVar_PositionInfo<-paste(subvarlist_forplot$Variant_Name_New_Cut, subvarlist_forplot$Variant.Details2, sep = ": ")
length(subvarlist_forplot$MainVar_PositionInfo)
x_order<-unique(subvarlist_forplot$MainVar_PositionInfo)
subvarlist_forplot$ordered_x<-factor(subvarlist_forplot$MainVar_PositionInfo, levels = x_order)
unique(subvarlist_forplot$MainVar_PositionInfo)

head(subvarlist_forplot)
unique(subvarlist_forplot$Variant.Details)
subvarlist_forplot[subvarlist_forplot$firstpos == "",]
subvarlist_forplot[subvarlist_forplot$firstpos == " ",] #this is pbs insertion, n=3
subvarlist_forplot[subvarlist_forplot$firstpos == ".1.",]


##organize and crunch main data to retain dates
submainlist$Date.of.collection.Prcsd<-lubridate::parse_date_time(submainlist$Date.of.collection,"ymd")
# #reorder by date 
submainlist<-submainlist[order(submainlist$Date.of.collection.Prcsd, decreasing = F),]
submainlist$rounded_week<-floor_date(submainlist$Date.of.collection.Prcsd, "week")

#get reordered date
submainlist$dateforgraph<-
  paste(month(submainlist$rounded_week),
        day(submainlist$rounded_week),
        year(submainlist$rounded_week), sep = "-" )

#get order then #force order  
countsperroundedweek_list<-unique(submainlist$dateforgraph)
countsperroundedweek_list_modified<-c(countsperroundedweek_list[1:2], "12-27-2020", "1-3-2021", "1-10-2021",
                                      countsperroundedweek_list[3:60], countsperroundedweek_list[61:115] )
submainlist$dateforgraph_factor<-factor(submainlist$dateforgraph, ordered=T, levels=countsperroundedweek_list_modified)

write.xlsx(submainlist, file = "samples_in_heatmap_Supp8v25data.xlsx")
write.xlsx(subvarlist_forplot, file = "variants_in_heatmap_Supp8v25.xlsx")
submainlist<-read.xlsx("samples_in_heatmap_Supp8v25data.xlsx")
subvarlist_forplot<-read.xlsx("variants_in_heatmap_Supp8v25.xlsx")
#need counts of all minor mutationss per week for both versions:
#counts per week
submainlist_perweek<-as.data.frame(submainlist %>% 
   group_by(dateforgraph_factor) %>% 
   dplyr::count(Variant.Details))
head(submainlist_perweek)
colnames(submainlist_perweek)<-c("dateforgraph_factor", "Variant.Details", "count_per_week")

#now merge my 2 lists
allcompiled_for8<-merge(submainlist_perweek, subvarlist_forplot, by.x = "Variant.Details", by.y = "Variant.Details")
head(allcompiled_for8)

#x axis will be dateforgraph_factor
#y axis will be ordered_x, fill will be count_per_week NOT n 

write.xlsx(allcompiled_for8, file = "SuppFig8_v25_fullheatmapdataset.xlsx")
allcompiled_for8<-read.xlsx("SuppFig8_v25_fullheatmapdataset.xlsx")

#switch date back to actual date 
allcompiled_for8$newDate<-as.Date(allcompiled_for8$dateforgraph_factor, format = "%m-%d-%Y")

allcompiled_for8<-allcompiled_for8[order(allcompiled_for8$ordervec, allcompiled_for8$firstpos, decreasing = FALSE),]

allcompiled_for8$ordered_x<-factor(allcompiled_for8$ordered_x, levels=unique(allcompiled_for8$ordered_x) )
allcompiled_for8$ordered_x<-forcats::fct_rev(allcompiled_for8$ordered_x)
head(allcompiled_for8,30)

#make new plot 
figSupp8black<-ggplot(allcompiled_for8, aes(x=newDate, y=ordered_x, fill = log10(count_per_week+1)))  +
    geom_tile() + ggtitle("All Subvariants ordered by Main Variant -> Position") + scale_fill_gradientn(colors=c(low = 'white',  high = 'black')) +
  xlab("Week") + ylab("MainVariant: SubVariant") + theme_bw() + scale_x_date(date_labels="%b %d %Y", date_breaks = "1 month") +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 5))

pdf("suppFig8_v25black.pdf", width = 10, height = 22)
figSupp8black
dev.off()

#vers 2 with gray - we used this version in the figure.
figSupp8gray<-ggplot(allcompiled_for8, aes(x=newDate, y=ordered_x, fill = log10(count_per_week+1)))  +
    geom_tile() + ggtitle("All Subvariants ordered by Main Variant -> Position") + scale_fill_gradientn(colors=c(low = 'lightgray',  high = 'black')) +
  xlab("Week") + ylab("MainVariant: SubVariant") + theme_bw() + scale_x_date(date_labels="%b %d %Y", date_breaks = "1 month") +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 5))

pdf("suppFig8_v25gray.pdf", width = 10, height = 22)
figSupp8gray
dev.off()
