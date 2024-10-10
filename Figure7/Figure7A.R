

#######final figure: fig 7 blue and purple heatmaps#
##make a new plot showing the VOC+ pos and the QS per voc

library(reshape2)

#load muttables
varfilewt<-read.csv("muttable_wt.csv")
varfilealpha<-read.csv("muttable_alpha.csv")
varfiledelta<-read.csv("muttable_delta.csv")
varfileomicron<-read.csv("muttable_omicron.csv")
varfileba2.75<-read.csv("muttable_omicronBA2.75.csv")
varfileba23<-read.csv("muttable_omicronBA23.csv")
varfileba45<-read.csv("muttable_omicronBA45.csv")
varfilexbb<-read.csv("muttable_XBB.csv")
varfilexbb.1.5<-read.csv("muttable_XBB.1.5.csv")
#some of these have multiple + in one 

#select everything after plus
varfilewt
varfilewt$variant<-"WT"
varfilewt$Freq<-NULL
#varfilewt$aachanges<-gsub("WT\\+","",varfilewt$WT.Variant)
varfilewt<-varfilewt %>% tidyr::separate_longer_delim(Var1, delim = c("+"))
varfilewt<-unique(varfilewt)
#get digits only
varfilewt$positionnum<-gsub("\\D+", "", varfilewt$Var1)
varfilewt

#alpha
varfilealpha
varfilealpha$variant<-gsub("\\+.*","",varfilealpha$Var1)
varfilealpha$Freq<-NULL; head(varfilealpha)
varfilealpha$aachanges<-gsub("B.1.1.7\\+","",varfilealpha$Var1)
varfilealpha<-varfilealpha %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfilealpha<-unique(varfilealpha)
#get digits only
varfilealpha$positionnum<-gsub("\\D+", "", varfilealpha$aachanges)
varfilealpha

##delta
varfiledelta
varfiledelta$variant<-gsub("\\+.*","",varfiledelta$Var1)
varfiledelta$Freq<-NULL; head(varfiledelta)
varfiledelta$aachanges<-gsub("B.1.617.2\\+","",varfiledelta$Var1)
varfiledelta<-varfiledelta %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfiledelta<-unique(varfiledelta)
#drop row with deletion
varfiledelta<-varfiledelta[varfiledelta$Var1!= "B.1.617.2+del23598-23600",]
#get digits only
varfiledelta$positionnum<-gsub("\\D+", "", varfiledelta$aachanges)
varfiledelta

###omicron 
varfileomicron
varfileomicron$variant<-gsub("\\+.*","",varfileomicron$Var1)
varfileomicron$Freq<-NULL; head(varfileomicron)
varfileomicron$aachanges<-gsub("B.1.1.529\\+","",varfileomicron$Var1)
varfileomicron<-varfileomicron %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfileomicron<-unique(varfileomicron)
#get digits only
varfileomicron$positionnum<-gsub("\\D+", "", varfileomicron$aachanges)
varfileomicron

###omicron bA 23
varfileba23
varfileba23$variant<-gsub("\\+.*","",varfileba23$Var1)
varfileba23$Freq<-NULL; head(varfileba23)
varfileba23$aachanges<-gsub("B.1.1.529 \\(BA.2\\/BA.3\\)\\+","",varfileba23$Var1)
varfileba23<-varfileba23 %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfileba23<-unique(varfileba23)
#get digits only
varfileba23$positionnum<-gsub("\\D+", "", varfileba23$aachanges)
varfileba23

###omicron ba 2.75.
varfileba2.75
varfileba2.75$variant<-gsub("\\+.*","",varfileba2.75$Var1)
varfileba2.75$Freq<-NULL; head(varfileba2.75)
varfileba2.75$aachanges<-gsub("B.1.1.529 \\(BA.2.75\\)\\+","",varfileba2.75$Var1)
varfileba2.75$aachanges<-gsub("B.1.1.529 \\(BA.2.75.2\\)\\+","",varfileba2.75$aachanges)
varfileba2.75[varfileba2.75$variant == "B.1.1.529 (BA.2.75.2)",]$variant<-"B.1.1.529 (BA.2.75)"
#had to do 2 because we have some as 2.75 and some as 2.75.2
varfileba2.75<-varfileba2.75 %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfileba2.75<-unique(varfileba2.75)
#dump row where var is B.1.1.529 (BA.2.75.2)
varfileba2.75<-varfileba2.75[varfileba2.75$Var1!="B.1.1.529 (BA.2.75.2)",]
varfileba2.75<-varfileba2.75[varfileba2.75$Var1!="B.1.1.529 (BA.2.75.2)+",]
#get digits only
varfileba2.75$positionnum<-gsub("\\D+", "", varfileba2.75$aachanges)
varfileba2.75
#one empty one because of the silent change


###omicron bA 45
varfileba45
varfileba45$variant<-gsub("\\+.*","",varfileba45$Var1)
varfileba45$Freq<-NULL; head(varfileba45)
varfileba45$aachanges<-gsub("B.1.1.529 \\(BA.4\\/BA.5\\)\\+","",varfileba45$Var1)
varfileba45<-varfileba45 %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfileba45<-unique(varfileba45)
varfileba45<-varfileba45[varfileba45$Var1!="B.1.1.529 (BA.4/BA.5)+",]
#get digits only
varfileba45$positionnum<-gsub("\\D+", "", varfileba45$aachanges)
varfileba45
#again a blank one due to the silent change 

###xbb
varfilexbb
varfilexbb$variant<-gsub("\\+.*","",varfilexbb$Var1)
varfilexbb$Freq<-NULL; head(varfilexbb)
varfilexbb$aachanges<-gsub("XBB\\+","",varfilexbb$Var1)
varfilexbb<-varfilexbb %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfilexbb<-unique(varfilexbb)
#get digits only
varfilexbb$positionnum<-gsub("\\D+", "", varfilexbb$aachanges)
varfilexbb<-varfilexbb[varfilexbb$positionnum != 15,]


###xbb1.5.
varfilexbb.1.5
varfilexbb.1.5$variant<-gsub("\\+.*","",varfilexbb.1.5$Var1)
varfilexbb.1.5$Freq<-NULL; head(varfilexbb.1.5)
varfilexbb.1.5$aachanges<-gsub("XBB.1.5\\+","",varfilexbb.1.5$Var1)
varfilexbb.1.5<-varfilexbb.1.5 %>% tidyr::separate_longer_delim(aachanges, delim = c("+"))
varfilexbb.1.5<-unique(varfilexbb.1.5)
#get digits only
varfilexbb.1.5$positionnum<-gsub("\\D+", "", varfilexbb.1.5$aachanges)
varfilexbb.1.5
varfilexbb.1.5<-varfilexbb.1.5[varfilexbb.1.5$positionnum != 15,]


colnames(varfilewt)<-c("AAChange","variant", "position")
varfilewt$subvariant<-varfilewt$AAChange
varfilewt<-varfilewt[,c("subvariant", "AAChange","variant", "position")]
colnames(varfilealpha)<-c("subvariant", "variant", "AAChange","position")
colnames(varfiledelta)<-c("subvariant", "variant", "AAChange", "position")
colnames(varfileomicron)<-c("subvariant", "variant", "AAChange", "position")
colnames(varfileba23)<-c("subvariant", "variant", "AAChange", "position")
colnames(varfileba2.75)<-c("subvariant","variant", "AAChange","position")
colnames(varfileba45)<-c("subvariant", "variant", "AAChange", "position")
colnames(varfilexbb)<-c("subvariant", "variant", "AAChange", "position")
colnames(varfilexbb.1.5)<-c("subvariant", "variant", "AAChange", "position")

vars_for_fig7_heatmap<-rbind(varfilewt, varfilealpha, varfiledelta, varfileomicron, varfileba23, varfileba2.75, varfileba45, varfilexbb, varfilexbb.1.5)
nrow(vars_for_fig7_heatmap)#411

#drop any that are RDRP changes
rdrp_approved<-read.table("rdrp_mutations.txt")
vars_for_fig7_heatmap<-vars_for_fig7_heatmap[!(vars_for_fig7_heatmap$AAChange %in% rdrp_approved$V1),]
nrow(vars_for_fig7_heatmap) #353



###
write.xlsx(vars_for_fig7_heatmap, "fig7_list.xlsx")
vars_for_fig7_heatmap<-read.xlsx("fig7_list.xlsx")

#drop empty ones and one with spbs deletion
vars_for_fig7_heatmap<-subset(vars_for_fig7_heatmap, vars_for_fig7_heatmap$position != "" )
vars_for_fig7_heatmap<-subset(vars_for_fig7_heatmap, vars_for_fig7_heatmap$position != "2359823600" )
nrow(vars_for_fig7_heatmap) #353
vars_for_fig7_heatmap$position<-as.numeric(vars_for_fig7_heatmap$position)
vars_for_fig7_heatmap$amplicon<-"none"
vars_for_fig7_heatmap


vars_for_fig7_heatmap[vars_for_fig7_heatmap$position < 505,]$amplicon<-"SRBM"
vars_for_fig7_heatmap[vars_for_fig7_heatmap$position > 671 ,]$amplicon<-"SPBS"

table(vars_for_fig7_heatmap$variant)
vars_for_fig7_heatmap_srbm<-subset(vars_for_fig7_heatmap, vars_for_fig7_heatmap$amplicon == "SRBM")
vars_for_fig7_heatmap_spbs<-subset(vars_for_fig7_heatmap, vars_for_fig7_heatmap$amplicon == "SPBS")

vars_for_fig7_heatmap_srbm<-vars_for_fig7_heatmap_srbm[,c("variant", "position", "amplicon")]
vars_for_fig7_heatmap_srbm<-unique(vars_for_fig7_heatmap_srbm)
vars_for_fig7_heatmap_spbs<-vars_for_fig7_heatmap_spbs[,c("variant", "position", "amplicon")]
vars_for_fig7_heatmap_spbs<-unique(vars_for_fig7_heatmap_spbs)

###need to load in the QS results...see folder Figure6CandSupplemental13-14 on github
srbmqs<-read.xlsx("fig6_andsupFig13-14/processed_srbd_AA_overall_table_ForFig7.xlsx")
spbsqs<-read.xlsx("fig6_andsupFig13-14/processed_spbs_aa_overall_table_forFig7.xlsx")

unique(srbmqs$sampleID)
unique(spbsqs$sampleID)
#convert lowercase to uppercase
srbmqs$sampleID<-toupper(srbmqs$sampleID)
length(unique(srbmqs$sampleID)) #228 samples in SRBM
#convert lowercase to uppercase
spbsqs$sampleID<-toupper(spbsqs$sampleID)
length(unique(spbsqs$sampleID)) #21 samples in SPBS
#is it the same samples in SRBM and SPBS?
unique(spbsqs$sampleID) %in% unique(srbmqs$sampleID)
unique(spbsqs[!(spbsqs$sampleID %in% srbmqs$sampleID),]$sampleID)
#these samples have an SPBS seq but no SRBM seq:
#"W62100962" "X61509058" "X61700640" "X61405493" "X61602521" "X61702023" "X61703692" "X61707656" "X61703575" "X61707650" "X61700492" "X61607558"
#so that is 12 of 21 that are not in SRBM list and 9 of 21 that are in SRBM list.
#228+12=240 different samples contributing


#need to know what proportion that is 
srbmqs$run_variant<-paste(srbmqs$run, srbmqs$variant, sep = "_")
srbmqs_samplerunonly<-srbmqs[,c("sampleID", "run_variant")]
srbmqs_samplerunonly<-unique(srbmqs_samplerunonly)
table(srbmqs_samplerunonly$run_variant)
#  139-141_Delta 139-141_Omicron           27_WT        31_Alpha           31_WT 
#             12              33              45              38             100 

spbsqs$run_variant<-paste(spbsqs$run, spbsqs$variant, sep = "_")
spbsqs_samplerunonly<-spbsqs[,c("sampleID", "run_variant")]
spbsqs_samplerunonly<-unique(spbsqs_samplerunonly)
table(spbsqs_samplerunonly$run_variant)
#  139-141_Delta 139-141_Omicron        31_Alpha           31_WT 
#              6              10               4               1

#27-WT input is 216/2 =  108
#31-WT input is 386/2 = 193
#31-Alpha input is 150/2 = 75
#139-141_Delta input is 386/2 = 193
#139-141_Omicron input is 1082 / 2 = 541 
#108+193+75+193+541=1110

#note that the srbd normal analysis data starts at 470 (Thr470), 471Glu, 472I, 473Y,  474 Gln, 
#so qs info starts at I so 472 so add 471 to it
##now make one that is QS data 
#and for spbs the position1 N is pos 679 so add 678 to the spbs positions

srbmqs<-srbmqs[,c("position", "change",  "variant")]
head(srbmqs)
srbmqs$position<-srbmqs$position+471
srbmqs<-unique(srbmqs)
table(srbmqs$variant)
srbmqs$variant<-paste(srbmqs$variant, "QS", sep = " ")
#need to drop ones that say no when the same pos is a yes 
srbmqs_yes<-subset(srbmqs, srbmqs$change == "yes")
srbmqs_no<-subset(srbmqs, srbmqs$change == "no")

srbmqs_yes$variant_change<-paste(srbmqs_yes$variant, srbmqs_yes$position, sep = "_")
srbmqs_no$variant_change<-paste(srbmqs_no$variant, srbmqs_no$position, sep = "_")
srbmqs_no<-subset(srbmqs_no, !(srbmqs_no$variant_change %in% srbmqs_yes$variant_change))

#spbs
spbsqs<-spbsqs[,c("position", "change",  "variant")]
head(spbsqs)
spbsqs$position<-spbsqs$position+678
spbsqs<-unique(spbsqs)
table(spbsqs$variant)
spbsqs$variant<-paste(spbsqs$variant, "QS", sep = " ")
#need to drop ones that say no when the same pos is a yes 
spbsqs_yes<-subset(spbsqs, spbsqs$change == "yes")
spbsqs_no<-subset(spbsqs, spbsqs$change == "no")

spbsqs_yes$variant_change<-paste(spbsqs_yes$variant, spbsqs_yes$position, sep = "_")
spbsqs_no$variant_change<-paste(spbsqs_no$variant, spbsqs_no$position, sep = "_")
spbsqs_no<-subset(spbsqs_no, !(spbsqs_no$variant_change %in% spbsqs_yes$variant_change))

#rebind
srbmqs<-rbind(srbmqs_yes, srbmqs_no); srbmqs$variant_change<-NULL; srbmqs
spbsqs<-rbind(spbsqs_yes, spbsqs_no); spbsqs$variant_change<-NULL; spbsqs

srbmqs$variant<-gsub(" QS","",srbmqs$variant)
spbsqs$variant<-gsub(" QS","",spbsqs$variant)

#combine voc and qs 
##adjust normal voc data
head(vars_for_fig7_heatmap_srbm,3)
head(vars_for_fig7_heatmap_spbs,3)
#make colnames match then rbind. rbind col opts too 
head(srbmqs,3)
head(spbsqs,3)

vars_for_fig7_heatmap_srbm$type<-"VOC+"
vars_for_fig7_heatmap_spbs$type<-"VOC+"
srbmqs$type<-"QS"
spbsqs$type<-"QS"

#adjust
vars_for_fig7_heatmap_srbm$amplicon<-NULL
vars_for_fig7_heatmap_srbm$change<-"yes"
head(vars_for_fig7_heatmap_srbm)
vars_for_fig7_heatmap_srbm<-vars_for_fig7_heatmap_srbm[,c("position", "change", "variant", "type" )]

allsrbm<-rbind(vars_for_fig7_heatmap_srbm, srbmqs)

table(allsrbm$variant)
allsrbm[allsrbm$variant == "B.1.1.529",]$variant <- "Omicron"
allsrbm[allsrbm$variant == "B.1.1.529 (BA.2.75)",]$variant <- "Omicron BA.2.75"
allsrbm[allsrbm$variant == "B.1.1.529 (BA.4/BA.5)",]$variant <- "Omicron BA.4/5"
allsrbm[allsrbm$variant == "B.1.1.529 (BA.2/BA.3)",]$variant <- "Omicron BA.2/3"
allsrbm[allsrbm$variant == "B.1.617.2",]$variant <- "Delta"
allsrbm[allsrbm$variant == "B.1.1.7",]$variant <- "Alpha"


allsrbm_wide <- spread(allsrbm, type, change)
allsrbm_wide
#fill NAs with no
allsrbm_wide[is.na(allsrbm_wide)] <- "no"

#adjust spbs
vars_for_fig7_heatmap_spbs$amplicon<-NULL
vars_for_fig7_heatmap_spbs$change<-"yes"
head(vars_for_fig7_heatmap_spbs)
vars_for_fig7_heatmap_spbs<-vars_for_fig7_heatmap_spbs[,c("position", "change", "variant", "type" )]

allspbs<-rbind(vars_for_fig7_heatmap_spbs, spbsqs)

table(allspbs$variant)
allspbs[allspbs$variant == "B.1.1.529",]$variant <- "Omicron"
allspbs[allspbs$variant == "B.1.1.529 (BA.2.75)",]$variant <- "Omicron BA.2.75"
allspbs[allspbs$variant == "B.1.1.529 (BA.4/BA.5)",]$variant <- "Omicron BA.4/5"
allspbs[allspbs$variant == "B.1.1.529 (BA.2/BA.3)",]$variant <- "Omicron BA.2/3"
allspbs[allspbs$variant == "B.1.617.2",]$variant <- "Delta"
allspbs[allspbs$variant == "B.1.1.7",]$variant <- "Alpha"
table(allspbs$variant)

allspbs_wide <- spread(allspbs, type, change)
allspbs_wide
#fill NAs with no
allspbs_wide[is.na(allspbs_wide)] <- "no"
allspbs_wide




####finalize calling####
allsrbm_wide$call<-"none"
allsrbm_wide[allsrbm_wide$QS == "yes" & allsrbm_wide$`VOC+` == "no",]$call<-"QS"
allsrbm_wide[allsrbm_wide$QS == "no" & allsrbm_wide$`VOC+` == "yes",]$call<-"VOC+"
allsrbm_wide[allsrbm_wide$QS == "yes" & allsrbm_wide$`VOC+` == "yes",]$call<-"VOC+ & QS"
table(allsrbm_wide$position)



allspbs_wide$call<-"none"
allspbs_wide[allspbs_wide$QS == "yes" & allspbs_wide$`VOC+` == "no",]$call<-"QS"
allspbs_wide[allspbs_wide$QS == "no" & allspbs_wide$`VOC+` == "yes",]$call<-"VOC+"
allspbs_wide[allspbs_wide$QS == "yes" & allspbs_wide$`VOC+` == "yes",]$call<-"VOC+ & QS"

###trim spare cols
head(allsrbm_wide)
allsrbm_wide<-allsrbm_wide[,c("position",  "variant", "call")]
allspbs_wide<-allspbs_wide[,c("position",  "variant",  "call")]


table(allsrbm_wide$variant)



unique(allsrbm_wide$variant)
allsrbm_wide$ordered_Variants_F<-factor(allsrbm_wide$variant, levels=c("XBB.1.5", "XBB", "Omicron BA.4/5", "Omicron BA.2.75",
                               "Omicron BA.2/3",  "Omicron", "Delta", "Alpha", "WT") )

unique(allspbs_wide$variant)
allspbs_wide$ordered_Variants_F<-factor(allspbs_wide$variant, levels=c("XBB.1.5", "XBB", "Omicron BA.4/5", "Omicron BA.2.75",
                               "Omicron BA.2/3",  "Omicron", "Delta", "Alpha", "WT") )


write.xlsx(allsrbm_wide, file = "fig7_allsrbm_heatmap.xlsx")
write.xlsx(allspbs_wide, file = "fig7_allspbs_heatmap.xlsx")

####finish this 
###make  vers with the voc plus and qs combo
library(tidyr)
allsrbm_wide<-read.xlsx("fig7_allsrbm_heatmap.xlsx")
head(allsrbm_wide)

table(allsrbm_wide$variant)
allsrbm_wide$ordered_Variants_F<-factor(allsrbm_wide$variant,levels=c("XBB.1.5", "XBB", "Omicron BA.4/5", "Omicron BA.2.75",
                               "Omicron BA.2/3",  "Omicron", "Delta", "Alpha", "WT"))
allsrbm_wide$position<-as.numeric(allsrbm_wide$position)

colopt5<-c("none"="white", "QS"="purple", "VOC+"="skyblue", "VOC+ & QS"="orchid1")

table(allsrbm_wide$ordered_Variants_F)

allsrbm_wide<-allsrbm_wide[allsrbm_wide$position>471,]
df2 <- allsrbm_wide %>% complete(ordered_Variants_F,position)

new_hm_vocs_qs_All_v2<-ggplot(df2, aes(x=position, y=ordered_Variants_F, fill = call)) +
    geom_tile(colour="black") + ggtitle("S-RBM Positions with VOC+ Changes and QS Hits") + 
  scale_fill_manual(values=colopt5, na.value = "white") +
  xlab("S-RBM Position") + ylab("VOC+/QS Hits") + theme_bw() + #scale_x_continuous()
  scale_x_continuous(limits = c(471,505), breaks = c(seq(472,504,1))) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 5,  hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 8))
new_hm_vocs_qs_All_v2

pdf("Fig7_SRBM_vocplusandqs.pdf", width = 6, height = 3)
new_hm_vocs_qs_All_v2
dev.off()


###repeat for spbs####
allspbs_wide<-read.xlsx("fig7_allspbs_heatmap.xlsx")
head(allspbs_wide)

table(allspbs_wide$variant)
allspbs_wide$ordered_Variants_F<-factor(allspbs_wide$variant,levels=c("XBB.1.5", "XBB", "Omicron BA.4/5", "Omicron BA.2.75",
                               "Omicron BA.2/3",  "Omicron", "Delta", "Alpha", "WT"))
allspbs_wide$position<-as.numeric(allspbs_wide$position)

colopt5<-c("none"="white", "QS"="purple", "VOC+"="skyblue", "VOC+ & QS"="orchid1")

table(allspbs_wide$ordered_Variants_F)
min(allspbs_wide$position)
max(allspbs_wide$position)

allspbs_wide<-allspbs_wide[allspbs_wide$position>678,]
df3 = allspbs_wide %>% complete(ordered_Variants_F,position)

new_hm_vocs_qs_spbs_v2<-ggplot(df3, aes(x=position, y=ordered_Variants_F, fill = call)) +
    geom_tile(colour="black") + ggtitle("S-PBS Positions with VOC+ Changes and QS Hits") + 
  scale_fill_manual(values=colopt5, na.value = "white") +
  xlab("S-PBS Position") + ylab("VOC+/QS Hits") + theme_bw() + #scale_x_continuous()
  scale_x_continuous(limits = c(678,699), breaks = c(seq(679,698,1))) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 5,  hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 8))
new_hm_vocs_qs_spbs_v2

pdf("Fig7_SPBS_vocplusandqs.pdf", width = 5, height = 3)
new_hm_vocs_qs_spbs_v2
dev.off()
