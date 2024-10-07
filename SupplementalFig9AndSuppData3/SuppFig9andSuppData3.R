#script for Supplemental Figure 9 and Supplemental Data 3. 
##Aggregating finalized variant calls into summary tables and generating summary barplots to compare the amplicons between main VOCs.

alldata<-read.xlsx("filtered_data_for_paper_v25.xlsx")
#fix dates
alldata$Date.of.collection<-as.Date(alldata$Date.of.collection, origin = "1899-12-30")
alldata$Date.of.collection.Prcsd<-lubridate::parse_date_time(alldata$Date.of.collection,"ymd")

#generate summary tables of the variant calls per main variant
table(alldata[alldata$Variant_Name_New == "WT",]$Variant.Details) 
table(alldata[alldata$Variant_Name_New == "WT+",]$Variant.Details) 
mut_table1<-table(alldata[alldata$Variant_Name_New == "WT+",]$Variant.Details) 
mut_table1<-mut_table1[order(mut_table1, decreasing = T)]
write.csv(mut_table1, "muttable_wt.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Alpha+",]$Variant.Details) 
mut_table2<-table(alldata[alldata$Variant_Name_New == "Alpha+",]$Variant.Details) 
mut_table2<-mut_table2[order(mut_table2, decreasing = T)]
write.csv(mut_table2, "muttable_alpha.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Delta+",]$Variant.Details) 
mut_table3<-table(alldata[alldata$Variant_Name_New == "Delta+",]$Variant.Details) 
mut_table3<-mut_table3[order(mut_table3, decreasing = T)]
write.csv(mut_table3, "muttable_delta.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Beta/Gamma+",]$Variant.Details) 
#only one here.
#P1/B.1.351+R683R 
#               9 

table(alldata[alldata$Variant_Name_New == "Omicron+",]$Variant.Details) 
mut_table4<-table(alldata[alldata$Variant_Name_New == "Omicron+",]$Variant.Details) 
mut_table4<-mut_table4[order(mut_table4, decreasing = T)]
write.csv(mut_table4, "muttable_omicron.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Omicron BA.2/3+",]$Variant.Details) 
mut_table5<-table(alldata[alldata$Variant_Name_New == "Omicron BA.2/3+",]$Variant.Details) 
mut_table5<-mut_table5[order(mut_table5, decreasing = T)]
write.csv(mut_table5, "muttable_omicronBA23.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Omicron BA.2.75+",]$Variant.Details) 
mut_table5.5<-table(alldata[alldata$Variant_Name_New == "Omicron BA.2.75+",]$Variant.Details) 
mut_table5.5<-mut_table5.5[order(mut_table5.5, decreasing = T)]
write.csv(mut_table5.5, "muttable_omicronBA2.75.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Omicron BA.4/5+",]$Variant.Details) 
mut_table6<-table(alldata[alldata$Variant_Name_New == "Omicron BA.4/5+",]$Variant.Details) 
mut_table6<-mut_table6[order(mut_table6, decreasing = T)]
write.csv(mut_table6, "muttable_omicronBA45.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "XBB+",]$Variant.Details) 
mut_table7<-table(alldata[alldata$Variant_Name_New == "XBB+",]$Variant.Details) 
mut_table7<-mut_table7[order(mut_table7, decreasing = T)]
write.csv(mut_table7, "muttable_XBB.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "XBB.1.5+",]$Variant.Details) 
mut_table8<-table(alldata[alldata$Variant_Name_New == "XBB.1.5+",]$Variant.Details) 
mut_table8<-mut_table8[order(mut_table8, decreasing = T)]
write.csv(mut_table8, "muttable_XBB.1.5.csv", quote = F, row.names = F)

table(alldata[alldata$Variant_Name_New == "Other",]$Variant.Details) 
#B.1.617.1 PBS insertion 
#            1             3 

mut_table9<-table(alldata[alldata$Variant_Name_New == "Other",]$Variant.Details) 
mut_table9<-mut_table9[order(mut_table9, decreasing = T)]
write.csv(mut_table9, "muttable_other.csv", quote = F, row.names = F)


table(alldata$Variant.Details) 
mut_tableAll<-table(alldata$Variant.Details)
mut_tableAll<-mut_tableAll[order(mut_tableAll, decreasing = T)]
write.csv(mut_tableAll, "muttable_alldata.csv", quote = F, row.names = F)
nrow(mut_tableAll)
#406
##we decided to compile the above muttables into an excel workbook and then ended up using them for further plots below

##Sup Data 3:
#use the same data to get counts and percents for Main variants and main variants+ for Supplemental Data 3
nrow(alldata) # 66165

##first tab - detailed variant counts
firsttab<-table(alldata$Variant.Details) 
firsttab<-firsttab[order(firsttab, decreasing = T)]

#second tab: main var percents (not subvars) so use Group
secondtab<-table(alldata$VariantGroup)
secondtab<-secondtab[order(secondtab, decreasing = T)]
secondtab<-secondtab/66165*100

#third tab: main var counts (not subvars) so use Group - same as 2nd but not percents
thirdtab<-table(alldata$VariantGroup)
thirdtab<-thirdtab[order(thirdtab, decreasing = T)]

##fourth tab: main vars percent but split into canonical and plus so use Variant Name New
fourthtab<-table(alldata$Variant_Name_New)
fourthtab<-fourthtab[order(fourthtab, decreasing = T)]
fourthtab<-fourthtab/66165*100

##fifth tab - same as 4th but using count
fifthtab<-table(alldata$Variant_Name_New)
fifthtab<-fifthtab[order(fifthtab, decreasing = T)]


list_of_datasets <- list("Detailed Variant Counts" = firsttab, "Main Variant Percents" = secondtab,
                         "Main Variant Counts" = thirdtab, "Main and Subvar Percents" = fourthtab, 
                         "Main and Subvar Counts" = fifthtab)
#write.xlsx(list_of_datasets, file = "supplData3_data-v25.xlsx")
#We also manually added canonical/non canonical labelling 

#load file - these are the manually compiled muttables from above.
varfile<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "All Variants")
varfilewt<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "WT Variants")
varfilealpha<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "Alpha Variants")
varfilebetagamma<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BetaGamma Variants")
varfiledelta<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "Delta Variants")
varfileomicron<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "Omicron Variants")
varfileba2.75<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BA2.75 Variants")
varfileba23<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BA23 Variants")
varfileba45<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BA45 Variants")
varfilexbb<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "XBB Variants")
varfilexbb.1.5<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "XBB.1.5 Variants")

varcounts<-as.data.frame(rbind(c("WT", nrow(varfilewt)), c("Alpha", nrow(varfilealpha)), c("Beta/Gamma", nrow(varfilebetagamma)), c("Delta", nrow(varfiledelta)), 
          c("Omicron", nrow(varfileomicron)), c("Omicron BA.2/3", nrow(varfileba23)), c("Omicron BA.2.75", nrow(varfileba2.75)),
          c("Omicron BA.4/5", nrow(varfileba45)), 
          c("XBB", nrow(varfilexbb)), c("XBB.1.5", nrow(varfilexbb.1.5))) )
colnames(varcounts)<-c("Variant", "NumUniqueVars")
varcounts$Variant<-factor(varcounts$Variant, levels=c("WT", "Alpha", "Beta/Gamma", "Delta", "Omicron", "Omicron BA.2/3", "Omicron BA.2.75", "Omicron BA.4/5",
                                                      "XBB", "XBB.1.5"))
varcounts$NumUniqueVars<-as.numeric(varcounts$NumUniqueVars)

##alpha	beta/gamma	delta	eta	mu	omicron	omicron 2.75.2	omicron 2/3	omicron 4/5	other	WT	XBB	XBB.1.5
#F8766D	#E18A00	#BE9C00	#8CAB00	#24B700	#00BE70	#00C1AB	#00BBDA	#00ACFC	#8B93FF	#D575FE	#F962DD	#FF65AC
# #alpha      b/g         delta.    eta         mu    omicron    2.75.    2/3       4/5       other     WT        XBB         XBB.1.5
#c("#F8766D" "#E18A00" "#BE9C00" "#8CAB00" "#24B700" "#00BE70" "#00C1AB" "#00BBDA" "#00ACFC" "#8B93FF" "#D575FE" "#F962DD" "#FF65AC")


colvec<-c("WT" = "#D575FE", "Alpha" = "#F8766D", "Beta/Gamma" = "#E18A00", "Delta" = "#BE9C00", "Omicron" = "#00BE70", "Omicron BA.2/3" = "#00BBDA", 
          "Omicron BA.2.75" = "#00C1AB", "Omicron BA.4/5" = "#00ACFC", "XBB" = "#F962DD", "XBB.1.5" = "#FF65AC")
colvec2<-c("WT" = "#D575FE", "B.1.1.7" = "#F8766D", "P1/B.1.351" = "#E18A00", "B.1.617.2" = "#BE9C00", "B.1.1.529" = "#00BE70", "B.1.1.529 (BA.2/BA.3)" = "#00BBDA", 
          "B.1.1.529 (BA.2.75)" = "#00C1AB", "B.1.1.529 (BA.4/BA.5)" = "#00ACFC", "XBB" = "#F962DD", "XBB.1.5" = "#FF65AC")


varbarplot<-ggplot(varcounts, aes(x=Variant, y=NumUniqueVars, fill = Variant))  +  geom_col(colour="black") + ggtitle("Unique Subvariants per Main Variant") + 
  scale_fill_manual(values=colvec) + xlab("Variant") + ylab("Count of Unique Subvariants") + theme_bw() + 
  theme(axis.title = element_text(size=10), legend.position = c("none"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),   axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
varbarplot

##next do panel of proportion of normal or voc+ per var
varfile<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "All Variants")
varfile
head(varfile)
#drop other, eta, mu
varfile<-varfile[!(varfile$variant %in% c("Other", "B.1.525", "PBS insertion", "B.1.617.1", "B.1.621" )) ,]
varfile$Group<-gsub("\\+.*", "", varfile$variant)
unique(varfile$Group)
#change 2.75.2 to be in same group as 2.75
varfile[varfile$Group == "B.1.1.529 (BA.2.75.2)",]$Group<-"B.1.1.529 (BA.2.75)"

varfile[!(varfile$Group %in% c("B.1.1.529","B.1.617.2","B.1.1.529 (BA.2/BA.3)", "B.1.1.7","B.1.1.529 (BA.4/BA.5)", "WT","P1/B.1.351","XBB.1.5" ,"B.1.1.529 (BA.2.75.2)",
                               "XBB",  "B.1.1.529 (BA.2.75)"  )),]$Group<-"WT"
varfile_pervar<-as.data.frame(varfile %>% 
    group_by(Group) %>% 
    dplyr::summarise(sumOfVar=sum(count)))

varfile<-merge(varfile, varfile_pervar, by = "Group")

varfile$VocPlus<-"VoC+"
varfile[varfile$Group == varfile$variant,]$VocPlus<-"VoC"

#do another round of group by to get just voc and voc+
varfile_pervar2<-as.data.frame(varfile %>% 
    group_by(Group, VocPlus) %>% 
    dplyr::summarise(sumOfVar2=sum(count)))
varfile_pervar3<-as.data.frame(varfile_pervar2 %>% 
    group_by(Group) %>% 
    dplyr::summarise(sumOfVar3=sum(sumOfVar2)))
varfileforgraph2<-merge(varfile_pervar2, varfile_pervar3, by= "Group")

varfileforgraph2$PercentOfVariant<-round(varfileforgraph2$sumOfVar2/varfileforgraph2$sumOfVar3 * 100, 2)

#remotes::install_github("coolbutuseless/ggpattern")
#library(ggpattern)
#library(ggplot2)

varfileforgraph2$Group<-factor(varfileforgraph2$Group, levels=c("WT", "B.1.1.7", "P1/B.1.351","B.1.617.2", "B.1.1.529", "B.1.1.529 (BA.2/BA.3)", "B.1.1.529 (BA.4/BA.5)", "B.1.1.529 (BA.2.75)",
                                                      "XBB", "XBB.1.5"))

varbarplot2<-ggplot(varfileforgraph2, aes(x=Group, y=PercentOfVariant, fill = Group, pattern = VocPlus))  +
    geom_col(position= position_stack()) + ggtitle("Proportion of VoC/Voc+ Per Main Variant") + 
  scale_fill_manual(values=colvec2) + scale_pattern_manual(values = c("VoC+" = "stripe", "VoC" = "none")) +
  geom_col_pattern(position= position_stack(), color = "black", pattern_fill = "black",  pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  xlab("Variant [solid = VoC, stripe = VoC+]") + ylab("Proportion") + theme_bw() + 
  theme(axis.title = element_text(size=10), legend.position = c("none"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text = element_text(size = 8),
      plot.title = element_text(hjust=0.5),  panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),   axis.text.y = element_text(hjust = 1, size = 10))


####third bars
varfile_mainvarcount<-unique(varfile_pervar3[,c("Group", "sumOfVar3")])
varfile_mainvarcount$Group<-factor(varfile_mainvarcount$Group, levels=c("WT", "B.1.1.7", "P1/B.1.351","B.1.617.2", "B.1.1.529", "B.1.1.529 (BA.2/BA.3)", "B.1.1.529 (BA.4/BA.5)", "B.1.1.529 (BA.2.75)",
                                                      "XBB", "XBB.1.5"))

##could also just add one starter bar plot which has just count of main vars incl voc+
varbarplot3<-ggplot(varfile_mainvarcount, aes(x=Group, y=sumOfVar3, fill = Group))  +
    geom_col(colour="black") + ggtitle("Count of Cases of Variants [as Log10]") + scale_fill_manual(values=colvec2) +
  xlab("Variant") + ylab("Log10 of Case Count") + theme_bw() + scale_y_log10(limits=c(1,100000), breaks=c(1,10,100,1000,10000,100000)) +
  theme(axis.title = element_text(size=10), legend.position = c("none"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text = element_text(size = 8),
      plot.title = element_text(hjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 10))

##fourth barplot showing number of unique AA spots altered in SRBD and PBS for the variant #
#so need table showing like alpha and which alpha spots had a +
varcounts
#load file
varfile<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "All Variants")
varfilewt<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "WT Variants")
varfilealpha<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "Alpha Variants")
varfilebetagamma<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BetaGamma Variants")
varfiledelta<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "Delta Variants")
varfileomicron<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "Omicron Variants")
varfileba2.75<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BA2.75 Variants")
varfileba23<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BA23 Variants")
varfileba45<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "BA45 Variants")
varfilexbb<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "XBB Variants")
varfilexbb.1.5<-read.xlsx("v25data_aggregatedvariants.xlsx", sheet = "XBB.1.5 Variants")


#Since rdrp and spbs amplicons have some overlaps we have a list of reference rdrp position to sort them out neatly for this plot
approvedrdrp<-read.table("rdrp_mutations.txt")
approvedrdrp
# 1     A685D
# 2     A685E
# 3     T686N
# 4     A685H
# 5     A685V
# 6     A699C
# 7  A702Stop
# 8     N705R
# 9     N705C
# 10 N713Stop
# 11    C697P
# 12    C697S
# 13 G712Stop
# 14    I696L
# 15    I715I
# 16    L707F
# 17    L708I
# 18    L708V
# 19    T701C
# 20    T701S
# 21    T710R
# 22    V700N
# 23    V700H
# 24    V700Y
# 25    V704K
# 26 V704Stop
# 27    V693F
# 28    D711W
#use this table to mark the rdrp changes
approvedrdrp<-approvedrdrp$V1


#write function to do crunch positions per amplicon for the variant 
splice_positions <- function(x,y){ #x = data table, y = name of variant 
  ##need count of pos that had a variant
  varfile<- x 
  varfile<-varfile$variant
  varfile<-c(strsplit(varfile, "\\+"))
  varfile_ind<-c()
  for(i in varfile){varfile_ind<-c(varfile_ind, i)}
  varfile_ind<-unique(varfile_ind) 
  #this is num of unique AA changes but we just want position numbers
  #convert to df and mark those which are RDRP approved 
  varfile_ind_df<-data.frame(varfile_ind, "amplicon"); colnames(varfile_ind_df)<-c("aachange", "amplicon")
  #get the position
  varfile_ind_df$position<-gsub("\\D+","",varfile_ind_df$aachange)
  #mark which Amplicon it belongs to
  varfile_ind_df[varfile_ind_df$position < 505,]$amplicon<-"SRBM"
  if(any(varfile_ind_df$position > 671)){varfile_ind_df[varfile_ind_df$position > 671,]$amplicon<-"SPBS"}
  #use reference list for RDRP
  if(any(varfile_ind_df$aachange %in% approvedrdrp)){varfile_ind_df[varfile_ind_df$aachange %in% approvedrdrp ,]$amplicon<-"RDRP"}
  varfile_ind_df_short<-varfile_ind_df[,2:3]
  #cut down to only amplicon and position columns then get unique only, then add variant type  
  varfile_ind_df_short<-unique(varfile_ind_df_short)
  varfile_ind_df_short$position<-as.numeric(varfile_ind_df_short$position)
  varfile_ind_df_short$variant<-y
  return(varfile_ind_df_short)
}

processed_wt<-splice_positions(varfilewt, "WT")

##for other variants need to preprocess by removing the main variant call from the front of the variant
varfilealpha$variant<-gsub("B.1.1.7\\+", "", varfilealpha$variant)
processed_alpha<-splice_positions(varfilealpha, "Alpha")

varfiledelta$variant<-gsub("B.1.617.2\\+", "", varfiledelta$variant)
#drop the deletion 
varfiledelta<-varfiledelta[varfiledelta$variant != "del23598-23600",]
processed_delta<-splice_positions(varfiledelta, "Delta")

varfilebetagamma$variant<-gsub("P1\\/B.1.351\\+", "", varfilebetagamma$variant)
#processed_betagamma<-splice_positions(varfilebetagamma, "Beta/Gamma")
#since this one is only R683R just manually adjust it 
processed_betagamma<-data.frame("SPBS", 683, "Beta/Gamma"); colnames(processed_betagamma)<-c("amplicon", "position", "variant")

varfileomicron$variant<-gsub("B.1.1.529\\+", "", varfileomicron$variant)
processed_omicron<-splice_positions(varfileomicron, "Omicron")

#for 2.75 drop the ones that are a ntd change only, and the ones that are 2.75.2
varfileba2.75$variant<-gsub("B.1.1.529 \\(BA.2.75\\)\\+", "", varfileba2.75$variant)
varfileba2.75$variant<-gsub("B.1.1.529 \\(BA.2.75.2\\)\\+", "", varfileba2.75$variant)
varfileba2.75$variant<-gsub("B.1.1.529 \\(BA.2.75.2\\)", "", varfileba2.75$variant)
varfileba2.75<-subset(varfileba2.75, varfileba2.75$variant != "")
processed_omicronba275<-splice_positions(varfileba2.75, "Omicron BA.2.75")

varfileba23$variant<-gsub("B.1.1.529 \\(BA.2\\/BA.3\\)\\+", "", varfileba23$variant)
processed_omicronba23<-splice_positions(varfileba23, "Omicron BA.2/3")

varfileba45$variant<-gsub("B.1.1.529 \\(BA.4\\/BA.5\\)\\+", "", varfileba45$variant)
varfileba45<-subset(varfileba45, varfileba45$variant  != "")
processed_omicronba45<-splice_positions(varfileba45, "Omicron BA.4/5")

varfilexbb$variant<-gsub("XBB\\+", "", varfilexbb$variant)
processed_xbb<-splice_positions(varfilexbb, "XBB")

varfilexbb.1.5$variant<-gsub("XBB.1.5\\+", "", varfilexbb.1.5$variant)
processed_xbb.1.5<-splice_positions(varfilexbb.1.5, "XBB.1.5")


##combine all of these things
allvarpositioncounts<-rbind(processed_wt, processed_alpha, processed_betagamma, processed_delta, processed_omicron, processed_omicronba275, processed_omicronba45, processed_xbb, processed_xbb.1.5)
#write.xlsx(allvarpositioncounts, file="uniquePositions_withVariants_perVariant_v25.xlsx")

#count per var and per amp
allvarpositioncounts_agg<-as.data.frame(allvarpositioncounts %>% 
    group_by(variant, amplicon) %>% 
      summarise(Positions = n_distinct(position)))

allvarpositioncounts_agg$Variant<-factor(allvarpositioncounts_agg$variant, levels = c("WT", "Alpha", "Beta/Gamma", "Delta", "Omicron", "Omicron BA.2/3", "Omicron BA.2.75", "Omicron BA.4/5", "XBB", "XBB.1.5"))
allvarpositioncounts_agg$Amplicon<-factor(allvarpositioncounts_agg$amplicon, levels=c("SRBM", "RDRP", "SPBS"))

##graph of mutated positions per variant######
varbarplot4<-ggplot(allvarpositioncounts_agg, aes(x=variant, y=Positions, fill = variant))  + facet_wrap(.~amplicon) +
    geom_col(colour="black", position='dodge2') + ggtitle("AA Positions per Amplicon with a Variant") + scale_fill_manual(values=colvec) +
  xlab("Variant") + ylab("Count") + theme_bw() +  scale_x_discrete(drop=FALSE) +
  theme(axis.title = element_text(size=10), legend.position = c("none"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text = element_text(size = 8),
      plot.title = element_text(hjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),    axis.text.y = element_text(hjust = 1, size = 10))
varbarplot4


pdf("barplotsofVarCounts_SupFig9_v25data.pdf", width = 5, height = 3)
varbarplot3
varbarplot
varbarplot2
varbarplot4
dev.off()

####end supp fig 9 bar plots#####
