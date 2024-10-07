


##Top part of panel A - the small tiles of antibody and ace2 info
#this is based on external data about ace2 affinity and antibody binding

data1<-read.xlsx("ab_with_srbd_heatmapdata.xlsx")
data1
data1<-subset(data1, data1$antibody %in% c("Absum", "HumanACE2", "VariantDefining"))

data1long<-melt(data1, id.vars=c("antibody"))
colnames(data1long)<-c("thing", "position", "hit")
data1long$position<-as.numeric(as.character(data1long$position))
thingorder<-c("Absum", "HumanACE2",  "VariantDefining")
data1long$thing<-factor(data1long$thing, levels = thingorder)
#drop variant defining, it was not needed in the final version
data1longsubset<-subset(data1long, data1long$thing != "VariantDefining")
data1longsubset<-subset(data1longsubset, data1longsubset$position < 505)

abhmmap2B<-ggplot(data1longsubset, aes(position, thing, fill = hit)) + geom_tile() +
  scale_fill_gradient(low="white",high = "navy", limits = c(0, 9), breaks = c(0,9),  labels = c(0,9)) +
    scale_x_continuous(limits = c(471,505), breaks = seq(472,504,1)) + ylab("") +
    theme(axis.title = element_text(size=5), legend.position = c("bottom"), legend.direction = ("horizontal"), plot.background = element_blank(), panel.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      axis.text.x = element_text(size = 6, hjust = 0.5), axis.text.y = element_text(size = 9))
abhmmap2B

pdf("2layer_boxes_v5_startsAt472ends504.pdf", width = 6, height = 2)
abhmmap2B
dev.off()

##The yellow/purple tiles of S-Pbs S1 and S2 are added in illustrator 


####lower part of 5A, the two red/white/black heatmaps of mutations at different AA positions####
re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v25.xlsx")
#fix dates
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")

###We first crunch the data and get AA per position per amplicon.

#call it data_for_minor variant plot, mvp 
data_for_mvp<-re_filtered_maindataset
data_for_mvp$Date.of.collection.Prcsd<-lubridate::parse_date_time(data_for_mvp$Date.of.collection,"ymd")
# #reorder by date 
data_for_mvp<-data_for_mvp[order(data_for_mvp$Date.of.collection.Prcsd, decreasing = F),]
week(data_for_mvp$Date.of.collection.Prcsd)
# https://lubridate.tidyverse.org/reference/round_date.html#rounding-up-date-objects-1
#try rounding down to week
data_for_mvp$rounded_week<-floor_date(data_for_mvp$Date.of.collection.Prcsd, "week")
#with this method every date is shifted to the sunday before.
#so a sat is shifted 6 days to the previous sunday

#get reordered date
data_for_mvp$dateforgraph<-
  paste(month(data_for_mvp$rounded_week),
        day(data_for_mvp$rounded_week),
        year(data_for_mvp$rounded_week), sep = "-" )
data_for_mvp$dateforgraph
#get order
countsperroundedweek_list<-unique(data_for_mvp$dateforgraph)
#get date list
data_for_mvp$Variant.Details

#need to look at all these separately 
table(data_for_mvp$Srbd_aa_var)
table(data_for_mvp$Srbd_v2_aa_var)

#fill in missing srbdv2 with srbd
data_for_mvp$srbddaas<-data_for_mvp$Srbd_v2_aa_var
data_for_mvp[is.na(data_for_mvp$srbddaas),]$srbddaas<-data_for_mvp[is.na(data_for_mvp$srbddaas),]$Srbd_aa_var
colnames(data_for_mvp)
#however, need to remove not-detected, Low-coverage warning labels from pipeline

data_for_mvp[is.na(data_for_mvp$Spbs_aa_var),]$sample
#assign everything with no value to be 'low'
data_for_mvp[is.na(data_for_mvp$srbddaas),]$srbddaas<-"low"
data_for_mvp[is.na(data_for_mvp$Spbs_aa_var),]$Spbs_aa_var<-"low"

#same for rdrp
data_for_mvp[is.na(data_for_mvp$RdRP_aa_var),]$sample # 
data_for_mvp[is.na(data_for_mvp$RdRP_aa_var),]$sample #~~2000 NA rdrps 
#mark as low
data_for_mvp[is.na(data_for_mvp$RdRP_aa_var),]$RdRP_aa_var<-"low"



#fixing by gsubbing - to "; " to match other rows 
data_for_mvp$srbddaas<-gsub("-", "; ", data_for_mvp$srbddaas)
data_for_mvp$Spbs_aa_var<-gsub("-", "; ", data_for_mvp$Spbs_aa_var)
#remove "; " from end if it is there
data_for_mvp$srbddaas<-gsub("; $", "", data_for_mvp$srbddaas)
data_for_mvp$Spbs_aa_var<-gsub("; $", "", data_for_mvp$Spbs_aa_var)

#we want to keep main var without + anno except for WT+ which has to still be called WT+ instead of other 
data_for_mvp$Variant_Name_New_All<-data_for_mvp$Variant_Name
data_for_mvp[data_for_mvp$Variant_Name_New == "WT+",]$Variant_Name_New_All<-"WT+"

table(data_for_mvp$Variant_Name)
table(data_for_mvp$Variant_Name_New)
table(data_for_mvp$Variant_Name_New_All)

data_for_mvp_srbd<-data_for_mvp[,c("sample","Variant.Details", "Date" ,"Variant_Name_New","Date.of.collection.Prcsd", 
                              "rounded_week", "dateforgraph", "srbddaas" )]
data_for_mvp_spbs<-data_for_mvp[,c("sample","Variant.Details", "Date" ,"Variant_Name_New","Date.of.collection.Prcsd", 
                              "rounded_week", "dateforgraph", "Spbs_aa_var" )]

#unnest
data_for_mvp_srbd_split<-data_for_mvp_srbd %>% separate_rows(srbddaas)
table(data_for_mvp_srbd_split$srbddaas)
#now table and sort the AA col
srbdtable<-table(data_for_mvp_srbd_split$srbddaas)
srbdtable<-srbdtable[order(srbdtable, decreasing = T)]
srbdtable
#need matchup of what indivual variant vs the variant name
srbdtable_x2<-table(data_for_mvp_srbd_split$srbddaas, data_for_mvp_srbd_split$Variant_Name_New)
#write.xlsx(srbdtable_x2, file="srbd_AAs_by_main_variant_v25.xlsx")

#repeat for spbs
data_for_mvp_spbs_split<-data_for_mvp_spbs %>% separate_rows(Spbs_aa_var)
#also an issue with V2 results where low was Low with L
table(data_for_mvp_spbs_split$Spbs_aa_var)
spbstable<-table(data_for_mvp_spbs_split$Spbs_aa_var)
spbstable<-spbstable[order(spbstable, decreasing = T)]
#need matchup of what indiv variant vs the variant name
spbstable_x2<-table(data_for_mvp_spbs_split$Spbs_aa_var, data_for_mvp_spbs_split$Variant_Name_New)
#write.xlsx(spbstable_x2, file="spbs_AAs_by_main_variant_v25.xlsx")

#Rdrp
#I added missing fields for RdRP_aa_var above already 
data_for_mvp$RdRP_aa_var
data_for_mvp_rdrp<-data_for_mvp[,c("sample","Variant.Details", "Date" ,"Variant_Name_New","Date.of.collection.Prcsd", 
                              "rounded_week", "dateforgraph", "RdRP_aa_var" )]
#fixing by gsubbing - to "; " to match other rows 
data_for_mvp_rdrp$RdRP_aa_var<-gsub("-", "; ", data_for_mvp_rdrp$RdRP_aa_var)
#remove "; " from end if RdRP_aa_var is there
data_for_mvp_rdrp$RdRP_aa_var<-gsub("; $", "", data_for_mvp_rdrp$RdRP_aa_var)
#unnest
data_for_mvp_rdrp_split<-data_for_mvp_rdrp %>% separate_rows(RdRP_aa_var)
data_for_mvp_rdrp_split[data_for_mvp_rdrp_split$RdRP_aa_var == "WTWT",]$RdRP_aa_var<-"WT" #typo from very early runs
rdrptable<-table(data_for_mvp_rdrp_split$RdRP_aa_var)
rdrptable<-rdrptable[order(rdrptable, decreasing = T)]
rdrptable
#need matchup of what indiv variant vs the variant name
rdrptable_x2<-table(data_for_mvp_rdrp_split$RdRP_aa_var, data_for_mvp_rdrp_split$Variant_Name_New)
#write.xlsx(rdrptable_x2, file="rdrp_AAs_by_main_variant_v25.xlsx")


### generate the copy of the data here
#fixing by gsubbing - to "; " to match other rows
data_for_mvp$srbddaas<-gsub("-", "; ", data_for_mvp$srbddaas)
data_for_mvp$Spbs_aa_var<-gsub("-", "; ", data_for_mvp$Spbs_aa_var)
#remove "; " from end if it is there
data_for_mvp$srbddaas<-gsub("; $", "", data_for_mvp$srbddaas)
data_for_mvp$Spbs_aa_var<-gsub("; $", "", data_for_mvp$Spbs_aa_var)
#fixing by gsubbing - to "; " to match other rows
data_for_mvp$RdRP_aa_var<-gsub("-", "; ", data_for_mvp$RdRP_aa_var)
#remove "; " from end if RdRP_aa_var is there
data_for_mvp$RdRP_aa_var<-gsub("; $", "", data_for_mvp$RdRP_aa_var)
table(data_for_mvp$Spbs_aa_var)

#unnest and save as three sheets
head(data_for_mvp)
data_for_mvp_srbd<-data_for_mvp[,c("sample","Variant.Details", "Date" ,"Variant_Name_New", "Date.of.collection.Prcsd",
                              "rounded_week", "dateforgraph", "srbddaas" )]
data_for_mvp_spbs<-data_for_mvp[,c("sample","Variant.Details","Date" ,"Variant_Name_New","Date.of.collection.Prcsd",
                              "rounded_week", "dateforgraph", "Spbs_aa_var" )]
data_for_mvp_rdrp<-data_for_mvp[,c("sample","Variant.Details", "Date" ,"Variant_Name_New","Date.of.collection.Prcsd",
                              "rounded_week", "dateforgraph", "RdRP_aa_var" )]

data_for_mvp_spbs_split<-data_for_mvp_spbs %>% separate_rows(Spbs_aa_var)
data_for_mvp_rdrp_split<-data_for_mvp_rdrp %>% separate_rows(RdRP_aa_var)
data_for_mvp_srbd_split<-data_for_mvp_srbd %>% separate_rows(srbddaas)

table(data_for_mvp_srbd_split$srbddaas)
table(data_for_mvp_spbs_split$Spbs_aa_var)
table(data_for_mvp_rdrp_split$RdRP_aa_var)

list_of_datasets <- list("SRBD" = data_for_mvp_srbd_split, "SPBS" = data_for_mvp_spbs_split, "RDRP" = data_for_mvp_rdrp_split)
write.xlsx(list_of_datasets, file = "minor_mutation_data_v25_includesWT+.xlsx")
#this file may occasionally differ from other aggregated tables of for example final variant calls
#since there are a few instances like in RDRP where long stretches of AA changes are detected but we suspect it is just a low quality amplicon result
#so those aren't reported to the final variant call of a sample


########fig 5 pink heatmaps ######
re_filtered_maindataset<-read.xlsx("filtered_data_for_paper_v25.xlsx")
#fix dates
re_filtered_maindataset$Date.of.collection<-as.Date(re_filtered_maindataset$Date.of.collection, origin = "1899-12-30")

library(RColorBrewer)

###Quick assessment of counts of vars
nrow(re_filtered_maindataset)
#66165
re_filtered_maindataset$Variant_Name_WTplus<-re_filtered_maindataset$Variant_Name_New
#changing this to include all WT counts
re_filtered_maindataset[re_filtered_maindataset$Variant_Name_New == "WT+",]$Variant_Name_WTplus<-"WT"
#should include all + here
re_filtered_maindataset$Variant_Name_WTplus<-gsub("\\+", "", re_filtered_maindataset$Variant_Name_WTplus)

table(re_filtered_maindataset$Variant_Name_WTplus)
#v25 data counts are used here 
#we exclude low count variants
totalcounts <- data.frame(variant = c("Alpha", "Beta/Gamma", "Delta", "Omicron", "Omicron BA.2/3", "Omicron BA.4/5", "WT", "XBB.1.5"),                    
                    counts = c(7807, 650, 9230, 27247,  8498, 7322,  4512, 568  )) #omicron was previously 29772 but now 27240 and BA45 was also slightly off 
totalcounts                    

#srbd
srbd_for_bp<-read.xlsx("minor_mutation_data_v25_includesWT+.xlsx", sheet = "SRBD")

srbd_for_bp$Variant_Name_New_All<-gsub("\\+","",srbd_for_bp$Variant_Name_New)
#fix wt+
srbd_for_bp[srbd_for_bp$Variant_Name_New == "WT+",]$Variant_Name_New_All<-"WT+"

srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$Variant_Name_New_All != "Other") #only 4 for other 
srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$Variant_Name_New_All != "Eta")
srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$Variant_Name_New_All != "Mu")
srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$Variant_Name_New_All != "Omicron BA.2.75")
srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$Variant_Name_New_All != "XBB")
srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$Variant_Name_New_All != "WT")

table(srbd_for_bp$Variant_Name_New_All)
#srbd_for_bp[srbd_for_bp$Variant_Name_New_All == "Beta",]$Variant_Name_New_All<-"Beta/Gamma"
srbd_for_bp<-subset(srbd_for_bp, srbd_for_bp$srbddaas != "WT")


#need to get the proper order by position so extract numbers from the mutations then add WT and low back on the end 
list_of_pos_srbd<-unique(srbd_for_bp$srbddaas)
list_of_pos_srbd_split1<-substring(list_of_pos_srbd, 0,3)
list_of_pos_srbd_split2<-substring(list_of_pos_srbd, 4,6)
list_of_pos_srbd_split2_num<-as.numeric(list_of_pos_srbd_split2)
list_of_pos_srbdsplit2_num<-list_of_pos_srbd_split2_num[!is.na(list_of_pos_srbd_split2_num)]
list_of_pos_srbd_split3<-substring(list_of_pos_srbd, 7,9)
ordered_srbd_<-seq(472, 504, 1)

srbd_for_bp$substring2<-substring(srbd_for_bp$srbddaas, 4,6)
#drop rows where it is WT which means this value is ""
srbd_for_bp<-subset(srbd_for_bp,(srbd_for_bp$substring2 != ""))
srbd_for_bp[srbd_for_bp$Variant_Name_New_All == "WT+",]
srbd_for_bp[srbd_for_bp$Variant_Name_New_All == "WT+",]$Variant_Name_New_All<-"WT"

#srbd_for_bp_per_main_var
srbd_for_bp_per_main_var<-as.data.frame(srbd_for_bp %>% 
   group_by(Variant_Name_New_All) %>% #could also drop + here to go back to normal variant name
   dplyr::count(substring2))

#now merge in total var nums
srbd_for_bp_per_main_var<-merge(srbd_for_bp_per_main_var, totalcounts, by.x = "Variant_Name_New_All", by.y = "variant", all =T)
srbd_for_bp_per_main_var$percent<-srbd_for_bp_per_main_var$n / srbd_for_bp_per_main_var$counts * 100

#w don't plot the 470 and 471 positions as they are in the primer.
srbd_for_bp_per_main_var<-subset(srbd_for_bp_per_main_var, !(srbd_for_bp_per_main_var$substring2 %in% c("470", "471")))
srbd_for_bp_per_main_var$substring2factor<-factor(srbd_for_bp_per_main_var$substring2, levels=c(ordered_srbd_))

table(srbd_for_bp_per_main_var$Variant_Name_New_All)

srbd_for_bp_per_main_var$Variant_Name_New_F<-factor(srbd_for_bp_per_main_var$Variant_Name_New_All, 
                                                     levels = c("WT", "Alpha", "Beta/Gamma","Delta","Omicron", "Omicron BA.2/3", "Omicron BA.4/5", "XBB.1.5"))

#set all canonical srbd location as NA so they will be gray on graph
srbd_for_bp_per_main_var_nc<-srbd_for_bp_per_main_var
srbd_for_bp_per_main_var_nc<-srbd_for_bp_per_main_var_nc[,c("percent", "substring2factor", "Variant_Name_New_F")]
#fill empty data
srbd_for_bp_per_main_var_nc<-as.data.frame(complete(srbd_for_bp_per_main_var_nc, Variant_Name_New_F, substring2factor))
#set these as 0
srbd_for_bp_per_main_var_nc[is.na(srbd_for_bp_per_main_var_nc$percent),]$percent<-0
#set canonical srbd spots to NA
srbd_for_bp_per_main_var_nc[is.na(srbd_for_bp_per_main_var_nc$percent),]
srbd_for_bp_per_main_var_nc[!(is.na(srbd_for_bp_per_main_var_nc$percent)) & srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "Alpha" & srbd_for_bp_per_main_var_nc$substring2factor == 501,]$percent<-NA
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "Beta/Gamma" & srbd_for_bp_per_main_var_nc$substring2factor %in% c(484, 501),]$percent<-NA
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "Delta" & srbd_for_bp_per_main_var_nc$substring2factor == 478,]$percent<-NA
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "Omicron" & 
                               srbd_for_bp_per_main_var_nc$substring2factor %in% c(477, 478, 484, 493, 496, 498, 501 ),]$percent<-NA
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "Omicron BA.2/3" & 
                               srbd_for_bp_per_main_var_nc$substring2factor %in% c(477, 478, 484, 493, 498, 501 ),]$percent<-NA
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "Omicron BA.4/5" & 
                               srbd_for_bp_per_main_var_nc$substring2factor %in% c(477, 478, 484, 486, 498, 501 ),]$percent<-NA
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$Variant_Name_New_F == "XBB.1.5" & 
                               srbd_for_bp_per_main_var_nc$substring2factor %in% c(477, 478, 484,486, 490, 498, 501 ),]$percent<-NA

#completed red
srbd_for_bp_per_main_var_nc$percent
hist(srbd_for_bp_per_main_var_nc$percent)
max(na.omit(srbd_for_bp_per_main_var_nc$percent)) #7.2

srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$percent > 5  & !(is.na(srbd_for_bp_per_main_var_nc$percent)) ,]
srbd_for_bp_per_main_var_nc[srbd_for_bp_per_main_var_nc$percent > 2.5  & !(is.na(srbd_for_bp_per_main_var_nc$percent)) ,]

#set the 2 values over 5 to just be 5 and label it as 5+
srbd_for_bp_per_main_var_nc_edited<-srbd_for_bp_per_main_var_nc
srbd_for_bp_per_main_var_nc_edited[srbd_for_bp_per_main_var_nc_edited$percent > 5  & !(is.na(srbd_for_bp_per_main_var_nc_edited$percent)) ,]$percent<-5


srbm_noncanonical<-ggplot(srbd_for_bp_per_main_var_nc_edited, aes(x=substring2factor, y=Variant_Name_New_F, fill = (percent) )) +
  ggtitle("Srbm Data: Variants, as percent of main variant total count (>5 set to 5)") + scale_x_discrete(drop=FALSE) +
  geom_tile(colour="black") + xlab("Position on Srbm") + ylab("Variant") + labs(fill="%") + theme_bw() +#
  scale_fill_gradientn(colours=c("white","tomato","red","brown","black"), breaks = c(0,2.5,5), limits = c(0,5), na.value="gray") + #
  theme(axis.title = element_text(size=12), legend.position = c("bottom"), panel.background=element_rect(fill="white", colour="white"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5, size = 10),
      panel.border = element_rect(colour = "black", fill=NA, size=1), legend.key.width= unit(1, 'cm'),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 12)  )
srbm_noncanonical

pdf("srbm_position_heatmap_normalizedtoVariantCases_v25_red_border_noncanonical.pdf", height = 3, width = 7)
srbm_noncanonical
dev.off()


####spbs####
spbs_for_bp<-read.xlsx("minor_mutation_data_v25_includesWT+.xlsx", sheet = "SPBS")

spbs_for_bp$Variant_Name_New_All<-gsub("\\+","",spbs_for_bp$Variant_Name_New)
#fix wt+
spbs_for_bp[spbs_for_bp$Variant_Name_New == "WT+",]$Variant_Name_New_All<-"WT+"
spbs_for_bp[spbs_for_bp$Variant_Name_New_All == "Beta",]$Variant_Name_New_All<-"Beta/Gamma"

head(spbs_for_bp)
spbs_for_bp<-subset(spbs_for_bp, spbs_for_bp$Variant_Name_New_All != "Other")
spbs_for_bp<-subset(spbs_for_bp, spbs_for_bp$Variant_Name_New_All != "Eta")
spbs_for_bp<-subset(spbs_for_bp, spbs_for_bp$Variant_Name_New_All != "Mu")
spbs_for_bp<-subset(spbs_for_bp, spbs_for_bp$Variant_Name_New_All != "Omicron BA.2.75")
spbs_for_bp<-subset(spbs_for_bp, spbs_for_bp$Variant_Name_New_All != "XBB")
spbs_for_bp<-subset(spbs_for_bp, spbs_for_bp$Spbs_aa_var != "WT")

#need to get the proper order by position so extract numbers from the mutations then add WT and low back on the end 
list_of_pos_spbs<-unique(spbs_for_bp$Spbs_aa_var)
list_of_pos_spbs<-subset(list_of_pos_spbs, list_of_pos_spbs!="low")
list_of_pos_spbs_split1<-substring(list_of_pos_spbs, 0,3)
list_of_pos_spbs_split2<-substring(list_of_pos_spbs, 4,6)
list_of_pos_spbs_split2_num<-as.numeric(list_of_pos_spbs_split2)
list_of_pos_spbs_split2_num<-list_of_pos_spbs_split2_num[!is.na(list_of_pos_spbs_split2_num)]
list_of_pos_spbs_split3<-substring(list_of_pos_spbs, 7,9)
ordered_spbs_<-seq(672, 699, 1)

spbs_for_bp$substring2<-substring(spbs_for_bp$Spbs_aa_var, 4,6)

#drop rows where it is WT which means this value is ""
spbs_for_bp<-subset(spbs_for_bp,(spbs_for_bp$substring2 != ""))
spbs_for_bp[spbs_for_bp$Variant_Name_New_All == "WT+",]
spbs_for_bp[spbs_for_bp$Variant_Name_New_All == "WT+",]$Variant_Name_New_All<-"WT"

spbs_for_bp_per_main_var<-as.data.frame(spbs_for_bp %>% 
   group_by(Variant_Name_New_All) %>%  #group_by(substring2) %>% 
   dplyr::count(substring2))

#now merge in total var nums
spbs_for_bp_per_main_var<-merge(spbs_for_bp_per_main_var, totalcounts, by.x = "Variant_Name_New_All", by.y = "variant", all =T)
spbs_for_bp_per_main_var$percent<-spbs_for_bp_per_main_var$n / spbs_for_bp_per_main_var$counts * 100

spbs_for_bp_per_main_var$substring2factor<-factor(spbs_for_bp_per_main_var$substring2, levels=c(ordered_spbs_))
spbs_for_bp_per_main_var$Variant_Name_New_F<-factor(spbs_for_bp_per_main_var$Variant_Name_New_All, 
                                                    levels = c("WT", "Alpha", "Beta/Gamma","Delta","Omicron", "Omicron BA.2/3", "Omicron BA.4/5", "XBB.1.5"))

####adjust complete cases, add NAs for canonical vars, fix scales if necessary
#set all canonical S26 location as NA so they will be gray on graph
spbs_for_bp_per_main_var_nc<-spbs_for_bp_per_main_var
spbs_for_bp_per_main_var_nc<-spbs_for_bp_per_main_var_nc[,c("percent", "substring2factor", "Variant_Name_New_F")]
#fill empty data
spbs_for_bp_per_main_var_nc<-as.data.frame(complete(spbs_for_bp_per_main_var_nc, Variant_Name_New_F, substring2factor))
#set these as 0
spbs_for_bp_per_main_var_nc[is.na(spbs_for_bp_per_main_var_nc$percent),]$percent<-0
#set canonical s26 spots to NA
spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "Alpha" & spbs_for_bp_per_main_var_nc$substring2factor == 681,]$percent<-NA
#spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "Beta/Gamma" & spbs_for_bp_per_main_var_nc$substring2factor %in% c(),]$percent<-NA
spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "Delta" & spbs_for_bp_per_main_var_nc$substring2factor == 681,]$percent<-NA
spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "Omicron" & 
                               spbs_for_bp_per_main_var_nc$substring2factor %in% c(679,681 ),]$percent<-NA
spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "Omicron BA.2/3" & 
                               spbs_for_bp_per_main_var_nc$substring2factor %in% c(679,681 ),]$percent<-NA
spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "Omicron BA.4/5" & 
                               spbs_for_bp_per_main_var_nc$substring2factor %in% c(679,681 ),]$percent<-NA
spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$Variant_Name_New_F == "XBB.1.5" & 
                               spbs_for_bp_per_main_var_nc$substring2factor %in% c(679,681 ),]$percent<-NA

spbs_for_bp_per_main_var_nc[spbs_for_bp_per_main_var_nc$percent > 5  & !(is.na(spbs_for_bp_per_main_var_nc$percent)) ,]
#highest value is just WT 681 at 12.44 %

spbs_for_bp_per_main_var_nc_capped<-spbs_for_bp_per_main_var_nc
spbs_for_bp_per_main_var_nc_capped[spbs_for_bp_per_main_var_nc_capped$percent > 5  & !(is.na(spbs_for_bp_per_main_var_nc_capped$percent)) ,]$percent<-5

spbssrbd_noncanonical_cap<-ggplot(spbs_for_bp_per_main_var_nc_capped, aes(x=substring2factor, y=Variant_Name_New_F, fill = (percent) )) +
  ggtitle("SPBS Data: Variants, as percent of main variant total count\n[>5 set to 5]") + scale_x_discrete(drop=FALSE) +
  geom_tile(colour="black") + xlab("Position on SPBS") + ylab("Variant") + labs(fill="%") + theme_bw() +#
  scale_fill_gradientn(colours=c("white","tomato","red","brown","black"), breaks = c(0,1,2,3,4,5), limits = c(0,5), na.value="gray") + #
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), panel.background=element_rect(fill="white", colour="white"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), legend.key.width= unit(1, 'cm'),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 12)  )
spbssrbd_noncanonical_cap

pdf("spbs_position_heatmap_normalizedtoVariantCases_v25_red_border_noncanonical_capped.pdf", height = 3.5, width = 7)
spbssrbd_noncanonical_cap
dev.off()

######End of 5A heatmaps ######
