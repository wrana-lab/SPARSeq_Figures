
####Next do Fig 7B Bloom data heatmaps.#####
##this analysis uses data from https://www.nature.com/articles/s41586-024-07636-1 by Bernadeta Dadonaite et al., Nature volume 631, pages617–626 (2024)

options(scipen=999)

###their 3 measurements are (fig 1 A)....
# BA.2 full-spike library with >9,000 mutations across spike
ba2spikeace2<-read.csv("bloompaper_spikeDMS_tables/bloompaper_Omicron_BA.2_spike_ACE2_binding_summary.csv")
head(ba2spikeace2)

# XBB.1.5 RBD-only saturated mutagenesis library
xbbrbddms<-read.csv("bloompaper_spikeDMS_tables/bloompaper_XBB.1.5_RBD_DMS_summary.csv")
head(xbbrbddms)

# XBB.1.5 full-spike library with >9,000 mutations across spike
xbbspikedms<-read.csv("bloompaper_spikeDMS_tables/bloompaper_XBB.1.5_spike_DMS_summary.csv")
head(xbbspikedms)

#####Start with SRBM######
####the by_main_variant files come from fig 5A script
fullsrbdsubvarlist<-read.xlsx("srbd_AAs_by_main_variant_v25.xlsx")

fullsrbdsubvarlist<-fullsrbdsubvarlist[fullsrbdsubvarlist$`Alpha+` > 0 | fullsrbdsubvarlist$`Beta/Gamma+` > 0 | fullsrbdsubvarlist$`Delta+` > 0 | fullsrbdsubvarlist$`Omicron+` > 0 | 
                                         fullsrbdsubvarlist$`XBB.1.5+` > 0 | fullsrbdsubvarlist$`XBB+` > 0 |
                                         fullsrbdsubvarlist$`Omicron.BA.2/3+` > 0 | fullsrbdsubvarlist$`Omicron.BA.2.75+` > 0 | fullsrbdsubvarlist$`Omicron.BA.4/5+` > 0 | fullsrbdsubvarlist$`WT+` > 0 ,]
nrow(fullsrbdsubvarlist) #97
head(fullsrbdsubvarlist)
fullsrbdsubvarlist_cut<-fullsrbdsubvarlist[,c(1,1)]
fullsrbdsubvarlist_cut<-subset(fullsrbdsubvarlist_cut, fullsrbdsubvarlist_cut$`.1` != "low")
fullsrbdsubvarlist_cut<-subset(fullsrbdsubvarlist_cut, fullsrbdsubvarlist_cut$`.1` != "WT")
nrow(fullsrbdsubvarlist_cut) #95
#cut up AAs to get original and mutant

fullsrbdsubvarlist_cut<-as.data.frame(fullsrbdsubvarlist_cut)
colnames(fullsrbdsubvarlist_cut)<-c("col1","col2")
#merge with qs list 
srbmqslist<-c("Gly496Ala", "Gly496Cys", "Gly496Arg", "Glu484Lys", "Ser477Asn", "Gly476Ser", "Pro479Leu", "Thr478Ile", "Asn501Tyr", "Ser477Ile", 
              "Ala475Val", "Val503Ile", "Ala475Ser", "Arg498Pro", "Thr500Ala")
srbmqslist<-data.frame(srbmqslist, srbmqslist)
colnames(srbmqslist)<-c("col1","col2")
fullsrbdsubvarlist_cut<-rbind(fullsrbdsubvarlist_cut,srbmqslist )
dim(fullsrbdsubvarlist_cut)
length(unique(fullsrbdsubvarlist_cut$col1))
fullsrbdsubvarlist_cut<-unique(fullsrbdsubvarlist_cut)

fullsrbdsubvarlist_cut$originalAA<-substr(fullsrbdsubvarlist_cut$col2, 1, 3)
fullsrbdsubvarlist_cut$position<-substr(fullsrbdsubvarlist_cut$col2, 4, 6)
fullsrbdsubvarlist_cut$mutantAA<-substr(fullsrbdsubvarlist_cut$col2, 7, 9)

#?a #a {seqinr}	R Documentation: Converts amino-acid three-letter code into the one-letter one
fullsrbdsubvarlist_cut$originalAA_1<-a(fullsrbdsubvarlist_cut$originalAA)
fullsrbdsubvarlist_cut$mutantAA_1<-a(fullsrbdsubvarlist_cut$mutantAA)

head(xbbspikedms)
xbbspikedms$sequential_site<-NULL
xbbspikedms$region<-NULL

#get combo column of wildtype, pos, mutant, to merge with our table
fullsrbdsubvarlist_cut$combocolumn<-paste0(fullsrbdsubvarlist_cut$position, fullsrbdsubvarlist_cut$mutantAA_1)
xbbspikedms$combocolumn<-paste0(xbbspikedms$site, xbbspikedms$mutant)

xbbspikedms_srbd<-subset(xbbspikedms, xbbspikedms$site < 505 & xbbspikedms$site > 469)
xbbspikedms_srbd$category<-"empty"
xbbspikedms_srbd[xbbspikedms_srbd$combocolumn %in% fullsrbdsubvarlist_cut$combocolumn,]$category<-"Normal"
xbbspikedms_srbd[!(xbbspikedms_srbd$combocolumn %in% fullsrbdsubvarlist_cut$combocolumn),]$category<-"BloomPaper"
xbbspikedms_srbd[xbbspikedms_srbd$human.sera.escape == 0 & xbbspikedms_srbd$spike.mediated.entry == 0 & xbbspikedms_srbd$ACE2.binding == 0 & xbbspikedms_srbd$category == "Normal",]$category<-"OurDataButZerosInBloomPaper"
table(xbbspikedms_srbd$category)

#find the missing part of subvar list cut where the combo isnt in the bloom table
fullsrbdsubvarlist_cut_toadd<-fullsrbdsubvarlist_cut[!(fullsrbdsubvarlist_cut$combocolumn %in% xbbspikedms_srbd$combocolumn),]
fullsrbdsubvarlist_cut_toadd$site<-fullsrbdsubvarlist_cut_toadd$position
fullsrbdsubvarlist_cut_toadd$mutant<-fullsrbdsubvarlist_cut_toadd$mutantAA_1
fullsrbdsubvarlist_cut_toadd$wildtype<-NA
fullsrbdsubvarlist_cut_toadd$human.sera.escape<-NA
fullsrbdsubvarlist_cut_toadd$spike.mediated.entry<-NA
fullsrbdsubvarlist_cut_toadd$ACE2.binding<-NA
fullsrbdsubvarlist_cut_toadd$category<-"OurPaper"
fullsrbdsubvarlist_cut_toadd<-fullsrbdsubvarlist_cut_toadd[,c("site", "wildtype", "mutant", "human.sera.escape", "spike.mediated.entry", "ACE2.binding", "combocolumn",   "category")]

#merge
srbdfinal<-rbind(xbbspikedms_srbd,fullsrbdsubvarlist_cut_toadd )
srbdfinal<-srbdfinal[srbdfinal$category!= "BloomPaper",]
srbdfinal$site<-as.numeric(srbdfinal$site)
srbdfinal$human.sera.escape<-as.numeric(srbdfinal$human.sera.escape)
srbdfinal$ACE2.binding<-as.numeric(srbdfinal$ACE2.binding)
srbdfinal$spike.mediated.entry<-as.numeric(srbdfinal$spike.mediated.entry)
write.xlsx(srbdfinal, file="srbd_bloom_results_fig7_cleanedvers.xlsx")

srbdfinal[srbdfinal$wildtype == srbdfinal$mutant & !(is.na(srbdfinal$wildtype)),]
srbdfinal$category2<-"Other"

#check
srbdfinal[srbdfinal$wildtype == srbdfinal$mutant & !(is.na(srbdfinal$wildtype)),]$category2<-"XBB.1.5 Canonical"
nrow(srbdfinal) #100
xbbfunctional<-subset(srbdfinal, srbdfinal$category2 != "Other")
xbbfunctional[xbbfunctional$combocolumn %in% fullsrbdsubvarlist_cut$combocolumn,]

#categorylist<-c("Normal" = NA, "OurDataButZerosInBloomPaper" = 8, "OurPaper" = 4) #"BloomPaper" = 3, 
category2list<-c("Other" = "gray", "XBB.1.5 Canonical" = "red") #"BloomPaper" = 3, 
write.xlsx(fullsrbdsubvarlist_cut, file="srbd_ourvocandqslist_results_fig7_cleanedvers_v2.xlsx")

#highlight those where the pos is higlighted (679, 681) and very high value
srbdfinal$serahighlight<-"none"
#srbdfinal[ (srbdfinal$site == 681 | srbdfinal$site == 679) & srbdfinal$human.sera.escape > 0 ,]$serahighlight<-"keyposition"
srbdfinal[!(is.na(srbdfinal$human.sera.escape)) & srbdfinal$human.sera.escape > 0 ,]$serahighlight<-"keyposition"
srbdfinal[srbdfinal$serahighlight != "none",]

srbdfinal$spikehighlight<-"none"
#srbdfinal[ (srbdfinal$site == 681 | srbdfinal$site == 679) & srbdfinal$spike.mediated.entry > 0 ,]$spikehighlight<-"keyposition"
srbdfinal[!(is.na(srbdfinal$spike.mediated.entry)) & srbdfinal$spike.mediated.entry > 0 ,]$spikehighlight<-"keyposition"
srbdfinal[srbdfinal$spikehighlight != "none",]

srbdfinal$ACE2.binding_highlight<-"none"
#srbdfinal[ (srbdfinal$site == 681 | srbdfinal$site == 679) & srbdfinal$ACE2.binding > 0 ,]$ACE2.binding_highlight<-"keyposition"
srbdfinal[!(is.na(srbdfinal$ACE2.binding)) & srbdfinal$ACE2.binding > 0 ,]$ACE2.binding_highlight<-"keyposition"
srbdfinal[srbdfinal$ACE2.binding_highlight != "none",]

highlightcols1 <- c("none" = "gray", keyposition = "red")

table(srbdfinal$serahighlight) #36 yes and 64 no
table(srbdfinal$spikehighlight) #28 yes and 72 no
table(srbdfinal$ACE2.binding_highlight) #13 yes and 13 no

#####repeat this for spbs########
fullspbssubvarlist<-read.xlsx("spbs_AAs_by_main_variant_v25.xlsx")

head(fullspbssubvarlist, 3)
#get only the + cols
fullspbssubvarlist<-fullspbssubvarlist[fullspbssubvarlist$`Alpha+` > 0 | fullspbssubvarlist$`Beta/Gamma+` > 0 | fullspbssubvarlist$`Delta+` > 0 | fullspbssubvarlist$`Omicron+` > 0 | 
                                         fullspbssubvarlist$`XBB.1.5+` > 0 | fullspbssubvarlist$`XBB+` > 0 |
                                         fullspbssubvarlist$`Omicron.BA.2/3+` > 0 | fullspbssubvarlist$`Omicron.BA.2.75+` > 0 | fullspbssubvarlist$`Omicron.BA.4/5+` > 0 | fullspbssubvarlist$`WT+` > 0 ,]
nrow(fullspbssubvarlist) #70
head(fullspbssubvarlist)
fullspbssubvarlist_cut<-fullspbssubvarlist[,c(1,1)]
fullspbssubvarlist_cut<-subset(fullspbssubvarlist_cut, fullspbssubvarlist_cut$`.1` != "low")
fullspbssubvarlist_cut<-subset(fullspbssubvarlist_cut, fullspbssubvarlist_cut$`.1` != "WT")
nrow(fullspbssubvarlist_cut) #68 but one is NA and gets dropped
#cut up AAs to get original and mutant

fullspbssubvarlist_cut<-as.data.frame(fullspbssubvarlist_cut)
colnames(fullspbssubvarlist_cut)<-c("col1","col2")
#merge with qs list 
spbsqslist<-c("Pro681His", "His681Pro", "Asn679Lys", "Arg681His", "Met696Lys", "Ser689Ile", "Lys679Asn")
spbsqslist<-data.frame(spbsqslist, spbsqslist)
colnames(spbsqslist)<-c("col1","col2")
fullspbssubvarlist_cut<-rbind(fullspbssubvarlist_cut,spbsqslist )


dim(fullspbssubvarlist_cut)
length(unique(fullspbssubvarlist_cut$col1))
fullspbssubvarlist_cut<-unique(fullspbssubvarlist_cut)
nrow(fullspbssubvarlist_cut)

fullspbssubvarlist_cut$originalAA<-substr(fullspbssubvarlist_cut$col2, 1, 3)
fullspbssubvarlist_cut$position<-substr(fullspbssubvarlist_cut$col2, 4, 6)
fullspbssubvarlist_cut$mutantAA<-substr(fullspbssubvarlist_cut$col2, 7, 9)
fullspbssubvarlist_cut<-subset(fullspbssubvarlist_cut, !(is.na(fullspbssubvarlist_cut$mutantAA)))
fullspbssubvarlist_cut<-subset(fullspbssubvarlist_cut, fullspbssubvarlist_cut$mutantAA!="NA")

fullspbssubvarlist_cut$originalAA_1<-a(fullspbssubvarlist_cut$originalAA)
fullspbssubvarlist_cut$mutantAA_1<-a(fullspbssubvarlist_cut$mutantAA)

nrow(fullspbssubvarlist_cut) #71
unique(fullspbssubvarlist_cut$combocolumn) #70 unique combo cols - 681H is twice 

xbbspikedms$sequential_site<-NULL
xbbspikedms$region<-NULL

#get combo column of wildtype, pos, mutant, to merge with our table
fullspbssubvarlist_cut$combocolumn<-paste0(fullspbssubvarlist_cut$position, fullspbssubvarlist_cut$mutantAA_1)
xbbspikedms$combocolumn<-paste0(xbbspikedms$site, xbbspikedms$mutant)

xbbspikedms_spbs<-subset(xbbspikedms, xbbspikedms$site < 700 & xbbspikedms$site > 670)
xbbspikedms_spbs$category<-"empty"
xbbspikedms_spbs[xbbspikedms_spbs$combocolumn %in% fullspbssubvarlist_cut$combocolumn,]$category<-"Normal"
xbbspikedms_spbs[!(xbbspikedms_spbs$combocolumn %in% fullspbssubvarlist_cut$combocolumn),]$category<-"BloomPaper"
xbbspikedms_spbs[xbbspikedms_spbs$human.sera.escape == 0 & xbbspikedms_spbs$spike.mediated.entry == 0 & xbbspikedms_spbs$ACE2.binding == 0 & xbbspikedms_spbs$category == "Normal",]$category<-"OurDataButZerosInBloomPaper"
table(xbbspikedms_spbs$category)

#find the missing part of subvar list cut where the combo isnt in the bloom table
fullspbssubvarlist_cut_toadd<-fullspbssubvarlist_cut[!(fullspbssubvarlist_cut$combocolumn %in% xbbspikedms_spbs$combocolumn),]
fullspbssubvarlist_cut_toadd$site<-fullspbssubvarlist_cut_toadd$position
fullspbssubvarlist_cut_toadd$mutant<-fullspbssubvarlist_cut_toadd$mutantAA_1
fullspbssubvarlist_cut_toadd$wildtype<-NA
fullspbssubvarlist_cut_toadd$human.sera.escape<-NA
fullspbssubvarlist_cut_toadd$spike.mediated.entry<-NA
fullspbssubvarlist_cut_toadd$ACE2.binding<-NA
fullspbssubvarlist_cut_toadd$category<-"OurPaper"
fullspbssubvarlist_cut_toadd<-fullspbssubvarlist_cut_toadd[,c("site", "wildtype", "mutant", "human.sera.escape", "spike.mediated.entry", "ACE2.binding", "combocolumn",   "category")]
head(xbbspikedms_spbs)

#merge
spbsfinal<-rbind(xbbspikedms_spbs,fullspbssubvarlist_cut_toadd )

spbsfinal<-spbsfinal[spbsfinal$category!= "BloomPaper",]
spbsfinal$site<-as.numeric(spbsfinal$site)
spbsfinal$human.sera.escape<-as.numeric(spbsfinal$human.sera.escape)
spbsfinal$ACE2.binding<-as.numeric(spbsfinal$ACE2.binding)
spbsfinal$spike.mediated.entry<-as.numeric(spbsfinal$spike.mediated.entry)
spbsfinal

write.xlsx(spbsfinal, file="spbs_bloom_results_fig7_cleanedvers.xlsx")
write.xlsx(fullspbssubvarlist_cut, file="spbs_ourlistofvocandqs_results_fig7_cleanedvers.xlsx")

spbsfinal$serahighlight<-"none"
#spbsfinal[ (spbsfinal$site == 681 | spbsfinal$site == 679) & spbsfinal$human.sera.escape > 0 ,]$serahighlight<-"keyposition"
spbsfinal[!(is.na(spbsfinal$human.sera.escape)) & spbsfinal$human.sera.escape > 0 ,]$serahighlight<-"keyposition"
spbsfinal[spbsfinal$serahighlight != "none",]

spbsfinal$spikehighlight<-"none"
#spbsfinal[ (spbsfinal$site == 681 | spbsfinal$site == 679) & spbsfinal$spike.mediated.entry > 0 ,]$spikehighlight<-"keyposition"
spbsfinal[!(is.na(spbsfinal$spike.mediated.entry)) & spbsfinal$spike.mediated.entry > 0 ,]$spikehighlight<-"keyposition"
spbsfinal[spbsfinal$spikehighlight != "none",]


spbsfinal$ACE2.binding_highlight<-"none"
#spbsfinal[ (spbsfinal$site == 681 | spbsfinal$site == 679) & spbsfinal$ACE2.binding > 0 ,]$ACE2.binding_highlight<-"keyposition"
spbsfinal[!(is.na(spbsfinal$ACE2.binding)) & spbsfinal$ACE2.binding > 0 ,]$ACE2.binding_highlight<-"keyposition"
spbsfinal[spbsfinal$ACE2.binding_highlight != "none",]

highlightcols1 <- c("none" = "gray", keyposition = "red")

table(spbsfinal$serahighlight)#17 yes 53 no
table(spbsfinal$spikehighlight) #13 yes 57 no
table(spbsfinal$ACE2.binding_highlight) #22 yes 48 no 

#use common min and max from srbd and spbs to design the scales so they are on the same scale system and we can use the same legends in the fig.
max(na.omit(srbdfinal$human.sera.escape)); min(na.omit(srbdfinal$human.sera.escape));
#0.3323, -1.25
max(na.omit(srbdfinal$spike.mediated.entry)); min(na.omit(srbdfinal$spike.mediated.entry));
#0.1395, -3.101
max(na.omit(srbdfinal$ACE2.binding)); min(na.omit(srbdfinal$ACE2.binding));
#0.9326, -1.881

###design heatmaps
#use common min and max from srbd and spbs to design the scales so they are on the same scale system and we can use the same legends in the fig.
max(c(na.omit(srbdfinal$human.sera.escape), na.omit(spbsfinal$human.sera.escape)))
#0.3323
min(c(na.omit(srbdfinal$human.sera.escape), na.omit(spbsfinal$human.sera.escape)))
#-1.25
max(c(na.omit(srbdfinal$spike.mediated.entry), na.omit(spbsfinal$spike.mediated.entry)))
#0.1395  
min(c(na.omit(srbdfinal$spike.mediated.entry), na.omit(spbsfinal$spike.mediated.entry)))
# -3.101
max(c(na.omit(srbdfinal$ACE2.binding), na.omit(spbsfinal$ACE2.binding)))
# 0.9326
min(c(na.omit(srbdfinal$ACE2.binding), na.omit(spbsfinal$ACE2.binding)))
# -1.881


#srbm
bloom_srbm_human.sera.escape<-ggplot(srbdfinal, aes(x=site, y=mutant, fill = human.sera.escape, colour = serahighlight)) + 
    geom_tile(aes(colour = serahighlight),size=0.5 ) + scale_colour_manual(values = highlightcols1) + ggtitle("human sera escape scores on our subvariant + QS sites") + 
  scale_fill_gradient2(low="blue", mid="white", high="darkred",  midpoint=0, limits=range(c(na.omit(srbdfinal$human.sera.escape), na.omit(spbsfinal$human.sera.escape))), na.value = "darkgray") +
  xlab("S-RBM Position") + ylab("Mutant AA") + theme_bw() + scale_x_continuous(limits = c(471,505), breaks = c(seq(472,504,1))) + 
  geom_tile(data = srbdfinal[srbdfinal$serahighlight == "keyposition",], aes(colour = serahighlight),size=0.5 ) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(), legend.direction = "vertical",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 12,  hjust = 1, angle = 45),
      axis.text.y = element_text(hjust = 1, size = 12))
bloom_srbm_human.sera.escape

bloom_srbm_spike.mediated.entry<-ggplot(srbdfinal, aes(x=site, y=mutant, fill = spike.mediated.entry, colour = spikehighlight)) + 
    geom_tile(aes(colour = spikehighlight),size=0.5 ) + scale_colour_manual(values = highlightcols1) + ggtitle("spike-mediated entry scores on our subvariant + QS sites") + 
  scale_fill_gradient2(low="blue", mid="white", high="darkred",  midpoint=0, limits=range(c(na.omit(srbdfinal$spike.mediated.entry), na.omit(spbsfinal$spike.mediated.entry))), na.value = "darkgray") +
  xlab("S-RBM Position") + ylab("Mutant AA") + theme_bw() + scale_x_continuous(limits = c(471,505), breaks = c(seq(472,504,1))) +
  geom_tile(data = srbdfinal[srbdfinal$spikehighlight == "keyposition",], aes(colour = spikehighlight),size=0.5 ) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),legend.direction = "vertical",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 12,  hjust = 1, angle = 45),
      axis.text.y = element_text(hjust = 1, size = 12))
bloom_srbm_spike.mediated.entry

bloom_srbm_ACE2.binding<-ggplot(srbdfinal, aes(x=site, y=mutant, fill = ACE2.binding, colour = ACE2.binding_highlight)) + 
    geom_tile(aes(colour = ACE2.binding_highlight),size=0.5 ) + scale_colour_manual(values = highlightcols1) + ggtitle("ACE2 binding scores on our subvariant + QS sites") +
  scale_fill_gradient2(low="blue", mid="white", high="darkred",  midpoint=0, limits=range(c(na.omit(srbdfinal$ACE2.binding), na.omit(spbsfinal$ACE2.binding))), na.value = "darkgray") +
  xlab("S-RBM Position") + ylab("Mutant AA") + theme_bw() + scale_x_continuous(limits = c(471,505), breaks = c(seq(472,504,1))) + 
  geom_tile(data = srbdfinal[srbdfinal$ACE2.binding_highlight == "keyposition",], aes(colour = ACE2.binding_highlight),size=0.5 ) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),legend.direction = "vertical",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 12,  hjust = 1, angle = 45),
      axis.text.y = element_text(hjust = 1, size = 12))
bloom_srbm_ACE2.binding

pdf("Fig7_Additional_SRBM_BloomDataAndOurData_MatchingScales.pdf", width = 5, height = 6)
bloom_srbm_human.sera.escape
bloom_srbm_spike.mediated.entry
bloom_srbm_ACE2.binding
dev.off()


#spbs
bloom_spbs_human.sera.escape<-ggplot(spbsfinal, aes(x=site, y=mutant, fill = human.sera.escape, colour = serahighlight)) + 
    geom_tile(aes(colour = serahighlight),size=0.5 ) + scale_color_manual(values=highlightcols1) + ggtitle("human sera escape scores on our subvariant + QS sites") +
  #intentionally using srbm range for scale
  scale_fill_gradient2(low="blue", mid="white", high="darkred",  midpoint=0, limits=range(c(na.omit(srbdfinal$human.sera.escape), na.omit(spbsfinal$human.sera.escape))), na.value = "darkgray") +
  xlab("S-PBS Position") + ylab("Mutant AA") + theme_bw() + scale_x_continuous(limits = c(678,699), breaks = c(seq(679,698,1))) + 
  geom_tile(data = spbsfinal[spbsfinal$serahighlight == "keyposition",], aes(colour = serahighlight),size=0.5 ) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(), legend.direction = "vertical",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 12,  hjust = 1, angle = 45),
      axis.text.y = element_text(hjust = 1, size = 12))
bloom_spbs_human.sera.escape

bloom_spbs_spike.mediated.entry<-ggplot(spbsfinal, aes(x=site, y=mutant, fill = spike.mediated.entry, colour = spikehighlight)) +
    geom_tile(aes(colour = spikehighlight),size=0.5 ) + scale_color_manual(values=highlightcols1) + ggtitle("spike-mediated entry scores on our subvariant + QS sites") +
  #intentionally using srbm range for scale
  scale_fill_gradient2(low="blue", mid="white", high="darkred",  midpoint=0, limits=range(c(na.omit(srbdfinal$spike.mediated.entry), na.omit(spbsfinal$spike.mediated.entry))), na.value = "darkgray") +
  xlab("S-PBS Position") + ylab("Mutant AA") + theme_bw() + scale_x_continuous(limits = c(678,699), breaks = c(seq(679,698,1))) + 
  geom_tile(data = spbsfinal[spbsfinal$spikehighlight == "keyposition",], aes(colour = spikehighlight),size=0.5 ) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),legend.direction = "vertical",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 12,  hjust = 1, angle = 45),
      axis.text.y = element_text(hjust = 1, size = 12))
bloom_spbs_spike.mediated.entry

bloom_spbs_ACE2.binding<-ggplot(spbsfinal, aes(x=site, y=mutant, fill = ACE2.binding, colour = ACE2.binding_highlight)) +
    geom_tile(aes(colour = ACE2.binding_highlight),size=0.5 ) + scale_color_manual(values=highlightcols1) + ggtitle("ACE2 binding scores on our subvariant + QS sites") + 
  #intentionally using srbm range for scale
  scale_fill_gradient2(low="blue", mid="white", high="darkred",  midpoint=0, limits=range(c(na.omit(srbdfinal$ACE2.binding), na.omit(spbsfinal$ACE2.binding))), na.value = "darkgray") +
  xlab("S-PBS Position") + ylab("Mutant AA") + theme_bw() + scale_x_continuous(limits = c(678,699), breaks = c(seq(679,698,1))) + 
  geom_tile(data = spbsfinal[spbsfinal$ACE2.binding_highlight == "keyposition",], aes(colour = ACE2.binding_highlight),size=0.5 ) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),legend.direction = "vertical",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(size = 12,  hjust = 1, angle = 45),
      axis.text.y = element_text(hjust = 1, size = 12))
bloom_spbs_ACE2.binding

pdf("Fig7_Additional_SPBS_BloomDataAndOurData_MatchingScales.pdf", width = 5, height = 6)
bloom_spbs_human.sera.escape
bloom_spbs_spike.mediated.entry
bloom_spbs_ACE2.binding
dev.off()
