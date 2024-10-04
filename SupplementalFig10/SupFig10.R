#Script to select random samples and plot some top and non top read percents for Supplemental Figure 10 
#this basically looked like same no matter what samples were selected 

re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v25.xlsx")
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")
head(re_filtered_maindataset)

####Get useful columns
#for V1 samples: "Srbd_total_read",	"Srbd_top1_read"
#for V2 samples: "Srbd_v2_total_read",	Srbd_v2_top1_read
table(re_filtered_maindataset$Variant_Name_New)
table(re_filtered_maindataset$VariantGroup)

#select from each main variant so
table(re_filtered_maindataset$Variant_Name_New)

sample_subset<-c()

for(i in c("Alpha", "Delta", "Omicron", "WT", "Alpha+", "Delta+", "Omicron+", "WT+", "Omicron BA.2.75",
           "Omicron BA.2.75+",   "Omicron BA.2/3",  "Omicron BA.2/3+",   "Omicron BA.4/5",  "Omicron BA.4/5+", "XBB.1.5", "XBB.1.5+")){
  ##NOT enough XBB or beta/gamma in our data to include those ones.
  set1<-subset(re_filtered_maindataset, re_filtered_maindataset$Variant_Name_New == i)
  if (i == "WT" | i == "WT+") {
      set1<-set1[!(is.na(set1$Srbd_total_read)),] 
  }
  if (nrow(set1) < 50) {
    sample_subset<-c(sample_subset, set1$sample) }

  else {miniset<-sample(set1$sample, 50, replace = FALSE)
    sample_subset<-c(sample_subset, miniset) }
}

setforgraph<-subset(re_filtered_maindataset, re_filtered_maindataset$sample %in% sample_subset)
table(setforgraph$Variant_Name_New)
#all 50 but XBB.1.5+ with only 26
setforgraph<-setforgraph[,c("sample", "Srbd_total_read",	"Srbd_top1_read", "Srbd_v2_total_read",	"Srbd_v2_top1_read", "Variant.Details", "Variant_Name_New") ]
class(setforgraph$Srbd_v2_total_read)
class(setforgraph$Srbd_v2_top1_read)
setforgraph$Srbd_v2_total_read<-as.numeric(setforgraph$Srbd_v2_total_read)
setforgraph$Srbd_v2_top1_read<-as.numeric(setforgraph$Srbd_v2_top1_read)
setforgraph$Srbd_v2_total_read<-round(setforgraph$Srbd_v2_total_read)

#organize percents
setforgraph$Srbd_top1_read_percent<-setforgraph$Srbd_top1_read / setforgraph$Srbd_total_read * 100
setforgraph$Srbd_other_read_percent<-100 - setforgraph$Srbd_top1_read_percent
setforgraph$Srbd_v2_top1_read_percent<-setforgraph$Srbd_v2_top1_read / setforgraph$Srbd_v2_total_read * 100
setforgraph$Srbd_v2_other_read_percent<-100 - setforgraph$Srbd_v2_top1_read_percent

#merge
setforgraph$top1_percent<-setforgraph$Srbd_top1_read_percent
setforgraph[is.na(setforgraph$top1_percent),]$top1_percent<-setforgraph[is.na(setforgraph$top1_percent),]$Srbd_v2_top1_read_percent
setforgraph$other_percent<-setforgraph$Srbd_other_read_percent
setforgraph[is.na(setforgraph$other_percent),]$other_percent<-setforgraph[is.na(setforgraph$other_percent),]$Srbd_v2_other_read_percent


##Make graph
setforgraphmini<-setforgraph[,c("top1_percent", "other_percent","Variant_Name_New")]
colnames(setforgraphmini)<-c("TopReadPercent", "OtherReadPercent", "Variant")
setforgraphmini<-melt(setforgraphmini, id.vars=c("Variant"))

#move WT to beginning
setforgraphmini$Variant_f = factor(setforgraphmini$Variant, levels=c("WT", "WT+", "Alpha","Alpha+","Delta","Delta+","Omicron","Omicron+",
                                                               "Omicron BA.2.75", "Omicron BA.2.75+", "Omicron BA.2/3", "Omicron BA.2/3+",
                                                               "Omicron BA.4/5",  "Omicron BA.4/5+", "XBB.1.5", "XBB.1.5+"))

#with top read % and non top read % on different graphs/pages
setforgraphmini_top<-subset(setforgraphmini, setforgraphmini$variable == "TopReadPercent")
setforgraphmini_nontop<-subset(setforgraphmini, setforgraphmini$variable == "OtherReadPercent")

setforgraphmini_new<-setforgraphmini
head(setforgraphmini_new)
setforgraphmini_new$VariantOnly<-gsub("\\+","",setforgraphmini_new$Variant)
setforgraphmini_new$PlusMinusOnly<-gsub("([a-zA-Z0-9])","",setforgraphmini_new$Variant)
setforgraphmini_new$PlusMinusOnly<-gsub("\\.","",setforgraphmini_new$PlusMinusOnly)


pdf("TopSRBMreadPlots_V25Data.pdf", width = 4, height = 5)
for(i in unique(setforgraphmini_new$VariantOnly)) {
  subset1<-subset(setforgraphmini_new, setforgraphmini_new$VariantOnly == i)
  subset1<-subset(subset1, subset1$variable == "TopReadPercent")
  print(ggplot(subset1, aes(x=Variant, y=value) ) + geom_boxplot() + ylim(0,100) +ggtitle("Srbm Top Read %\n[50 samples each except XBB.1.5+]") +
    xlab("Variant") + ylab("Percent") +   theme(axis.title = element_text(size=12), legend.position = c("none"), plot.background = element_blank(),
        plot.title = element_text(hjust=0.5, size = 8), strip.text.x = element_text(size=6), panel.border = element_rect(colour = "black", fill=NA, size=1), 
        strip.background = element_blank(), axis.text.x = element_text(size = 10, angle = 90, hjust = 1), axis.text.y = element_text(hjust = 1, size = 10)))
}
dev.off()
