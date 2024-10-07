#Script for Fig 5B - top subvariants heatmap and bar plot.

#Need to first make aggregated data of 'subvariants', so the +XYZ calls we assigned to samples.

##import v25 data
alldatareload<-read.xlsx("filtered_data_for_paper_v25.xlsx")
alldatareload$Date.of.collection<-as.Date(alldatareload$Date.of.collection, origin = "1899-12-30") #

###get count of all subvars and subvars per main var###
allsubvar<-table(alldatareload$Variant.Details)
allsubvar<-allsubvar[order(allsubvar, decreasing = T)]
write.xlsx(allsubvar, file="AllSubvariants_v25.xlsx")

##then break down by main var
#so group by variant name new and then variant details within that
alldatareload$Variant_Name_New_Cut<-gsub("\\+", "", alldatareload$Variant_Name_New)

data_for_bp<-as.data.frame(alldatareload %>%
   group_by(Variant_Name_New_Cut) %>%
   dplyr::count(Variant.Details))
#organized subvar data
write.xlsx(data_for_bp, "organized_subvariants_v25.xlsx")


###Begin bar plot of residue categories####
##Get top subvariants
subvarlist<-read.xlsx("organized_subvariants_v25.xlsx")

#need subset without the samples that are just canonical main VOCs
mainvars<-c("B.1.1.7", "B.1.1.529", "B.1.1.529 (BA.2/BA.3)", "B.1.1.529 (BA.4/BA.5)", "WT", "XBB", 
            "XBB.1.5", "B.1.1.525", "B.1.621", "B.1.617.2", "P1/B.1.351",
            "B.1.1.529 (BA.2.75)", "B.1.1.529 (BA.2.75.2)")

subvarlist<-subset(subvarlist, !(subvarlist$Variant.Details %in% mainvars))
#We want to select the top, most frequent variants in our data 
subvarlist_cut<-subset(subvarlist, subvarlist$n >= 40)
nrow(subvarlist_cut) #15 variants
subvarlist_cut$Variant.Details
write.xlsx(subvarlist_cut, "top14_subvariants_v25data.xlsx")


##import V25 data
alldatareload<-read.xlsx("filtered_data_for_paper_v25.xlsx")
alldatareload$Date.of.collection<-as.Date(alldatareload$Date.of.collection, origin = "1899-12-30") #

subvarlist<-read.xlsx("top14_subvariants_v25data.xlsx")
head(subvarlist)
#we picked these because they're the subvariants with the most counts ( > 40)

subvarlist<-c( "E484K", "E484K+P681H", "N501Y", "N501Y+P681R", "P681R", "P681H", "B.1.1.529 (BA.2.75)+F490S", "B.1.1.529+P681R", "B.1.1.529 (BA.4/BA.5)+A672G",
              "A684V", "S691S", "B.1.617.2+R682R", "B.1.617.2+A688V", "B.1.1.7+F490L") 

subsetdata<-subset(alldatareload, alldatareload$Variant.Details %in% subvarlist)
table(subsetdata$Variant.Details)
nrow(subsetdata)


##organize and crunch main data to retain dates
subsetdata$Date.of.collection.Prcsd<-lubridate::parse_date_time(subsetdata$Date.of.collection,"ymd")
# #reorder by date 
subsetdata<-subsetdata[order(subsetdata$Date.of.collection.Prcsd, decreasing = F),]
colnames(subsetdata)
subsetdata<-subsetdata[,c("sample","Variant.Details", "Variant_Name_New", "VariantGroup", "Date.of.collection", "Date.of.collection.Prcsd")]
#for each subvar need date of 1st one, so that we can calculate the number of days since 1st case for each subsequent case.

#initialize table and run loop
outtable<-data.frame()

for(i in unique(subsetdata$Variant.Details)){
  sub_table<-subset(subsetdata, subsetdata$Variant.Details == i)
  sub_table<-sub_table[c(1),]
  sub_table$CaseNum<-c(1)

  outtable<-rbind(outtable,sub_table)
}
outtable
class(outtable$Date.of.collection.Prcsd)
outtable<-outtable[,c("Variant.Details", "Date.of.collection.Prcsd")]
colnames(outtable)<-c("Variant.Details", "Date.of.first.case")

#merge first days 
subsetdata<-merge(subsetdata, outtable, by = "Variant.Details")
head(subsetdata)

#convert dates
subsetdata$NumDaysSince1stCase<-subsetdata$Date.of.collection.Prcsd - subsetdata$Date.of.first.case
subsetdata$NumDaysSince1stCase<-as.numeric(gsub(" days", "", subsetdata$NumDaysSince1stCase))
subsetdata$NumDaysSince1stCase<-subsetdata$NumDaysSince1stCase/86400
head(subsetdata)
subsetdata<-subsetdata[,c( "Variant_Name_New", "Variant.Details", "VariantGroup", "Date.of.collection.Prcsd", "Date.of.first.case", "NumDaysSince1stCase")]
#write.xlsx(subsetdata, file = "Top14SubvariantHeatmap_daysSinceFirstCaseOfSubvariant_v25data.xlsx")
subsetdata<-read.xlsx("Top14SubvariantHeatmap_daysSinceFirstCaseOfSubvariant_v25data.xlsx")

subsetdata$Date.of.collection.Prcsd<-NULL
subsetdata$Date.of.first.case<-NULL

submainlist_agg<-data.frame(subsetdata %>% group_by(Variant_Name_New, Variant.Details, NumDaysSince1stCase) %>% summarise(n = n()))
unique(submainlist_agg$Variant.Details)
# [1] "B.1.1.7+F490L"               "B.1.617.2+A688V"             "B.1.617.2+R682R"             "B.1.1.529 (BA.2.75)+F490S"   "B.1.1.529 (BA.4/BA.5)+A672G" "B.1.1.529+P681R"            
# [7] "A684V"                       "E484K"                       "E484K+P681H"                 "N501Y"                       "N501Y+P681R"                 "P681H"                      
#[13] "P681R"                       "S691S"   


#write.xlsx(submainlist_agg, "Top14SubvariantHeatmap_forMADcalc_v25data.xlsx")
submainlist_agg<-read.xlsx("Top14SubvariantHeatmap_forMADcalc_v25data.xlsx")

#we decided to split these into short and long lasting variants
subvar_order_short<-c( "B.1.1.529 (BA.4/BA.5)+A672G", "B.1.1.7+F490L",  "B.1.617.2+A688V", "B.1.617.2+R682R",  "S691S",  "A684V") #group 1
subvar_order_long<-c("B.1.1.529 (BA.2.75)+F490S",  "B.1.1.529+P681R", "P681H", "P681R",   "N501Y+P681R",  "N501Y", "E484K+P681H",
                 "E484K")
subvar_order<-c(subvar_order_short, subvar_order_long)

submainlist_short<-subset(subsetdata, subsetdata$Variant.Details %in% subvar_order_short)
submainlist_long<-subset(subsetdata, subsetdata$Variant.Details %in% subvar_order_long)

submainlist_agg$Variant.Details_F<-factor(submainlist_agg$Variant.Details, levels = subvar_order)
min(submainlist_agg$NumDaysSince1stCase) #0

#load functionality table
#need new table that is median days since first case
submainlist_agg %>%
  group_by(Variant.Details) %>%
  summarise(medianN = median(NumDaysSince1stCase))
#saved this manually and edited it to assign categories after discussion with the team
#then reloaded
forbarplot5b<-read.xlsx("MedianDaysSince1stCase_categorized_v25.xlsx")
forbarplot5b_other<-subset(forbarplot5b, forbarplot5b$Category == "Other residue(s)")
forbarplot5b_notother<-subset(forbarplot5b, forbarplot5b$Category != "Other residue(s)")

forbarplot5b$Variant.Details_F<-factor(forbarplot5b$Variant.Details, levels=c("B.1.1.529 (BA.4/BA.5)+A672G", "B.1.1.7+F490L", "B.1.617.2+A688V", "B.1.617.2+R682R", "S691S", 
                                                                              "A684V", "B.1.1.529 (BA.2.75)+F490S","B.1.1.529+P681R", "P681H", "P681R", "N501Y+P681R", "N501Y",
                                                                              "E484K+P681H", "E484K"))
forbarplot5b$Category_F<-factor(forbarplot5b$Category, levels = c("Variant-defining residue(s)", "Other residue(s)"))

forbarplot5b_plot<-ggplot(forbarplot5b, aes(x=medianN, y = Variant.Details_F, fill = Category_F)) + geom_col() + 
  ggtitle("5B median bars") +   theme_bw() + 
  theme(axis.title = element_text(size=0), legend.position = c("right"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 8,  hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 8))
forbarplot5b_plot
pdf("forbarplot5b_plot_median_v25.pdf", width=5,height=4)
forbarplot5b_plot
dev.off()

#Median of median values per category
median(forbarplot5b[forbarplot5b$Category == "Other residue(s)",]$medianN) #64.5
median(forbarplot5b[forbarplot5b$Category == "Variant-defining residue(s)",]$medianN) #29
median(forbarplot5b$medianN) #46
#We mention the MAD in the text,

###end of bar plot of residue categories

###begin heatmap of days since first case####
#check that each one actually has a value at 0
submainlist_agg[submainlist_agg$NumDaysSince1stCase<1,]
#write.xlsx(submainlist_agg, file = "data_for_heatmap5b_v25data.xlsx")
submainlist_agg<-read.xlsx("data_for_heatmap5b_v25data.xlsx")

keysubvar_plot_all<-ggplot(submainlist_agg, aes(x=NumDaysSince1stCase, y=Variant.Details_F, fill = n)) +
    geom_tile() + ggtitle("Key Subvariants: Days Since 1st Case") + scale_fill_gradientn(colors=c(low = 'thistle2',  high = 'black')) +
  xlab("Days Since 1st Case of Subvariant") + ylab("Subvariant") + theme_bw() + 
  guides(fill = guide_colourbar(title = "Number of cases")) + scale_x_continuous(limits = c(-2,200), breaks = c(0,50,100,150,200)) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 12,  hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 8))
keysubvar_plot_all


pdf("heatmap5B_timebased_topsubvariants_datav25.pdf", width = 8, height = 4)
keysubvar_plot_all
dev.off()

###End of 5b####
