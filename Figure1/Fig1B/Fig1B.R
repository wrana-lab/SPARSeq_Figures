##import V25 data
alldatareload<-read.xlsx("filtered_data_for_paper_v25.xlsx")
alldatareload$Date.of.collection<-as.Date(alldatareload$Date.of.collection, origin = "1899-12-30") 

# #####Rfinal version of pilot heatmap######
# ####Use the Nov 29 to Jan 28 2021 pilot data###
# #now with extra preclinical filter it is 1199 samples

#load list of approved samples from pilot time period
pilotsamples<-read.xlsx("nov2020tojan2021_pilotset_v3_sampleIDandRunIDs.xlsx")
head(pilotsamples)

pilotset<-subset(alldatareload, alldatareload$sample %in% pilotsamples$sample)
nrow(pilotset) #1198

firstfewmonths<-pilotset
#write.xlsx(firstfewmonths, file = "v25data_pilotset_forheatmap.xlsx")
#firstfewmonths<-read.xlsx("v25data_pilotset_forheatmap.xlsx")

#order by RBD zscore (rank order)##
firstfewmonths$testSrbd<- log10( (firstfewmonths$Srbd+1) / (firstfewmonths$ACTB+1) )
firstfewmonths<-firstfewmonths[order(firstfewmonths$testSrbd),]
orderedbyzscore<-as.factor(firstfewmonths$sample)
firstfewmonths_hm3<-firstfewmonths[,c("sample", "ACTB",	"Rdrp",	"Spoly",	"Srbd")]
pilotnegs<-read.xlsx("pilotset_negatives.xlsx", sheet = "forR")
pilotnegs<-pilotnegs[,c("sample",	"ACTB",	"Rdrp",	"Spoly",	"Srbd")]
firstfewmonths_hm3<-rbind(firstfewmonths_hm3, pilotnegs)
#firstfewmonths_hm3$sample<-rownames(firstfewmonths_hm3)
rownames(firstfewmonths_hm3)<-firstfewmonths_hm3$sample
firstfewmonths_hm3$sample<-NULL

#get actb ratio
# log10[(viral gene counts+1)/(actb counts+1)]
firstfewmonths_hm3$Rdrp<- log10( (firstfewmonths_hm3$Rdrp+1) / (firstfewmonths_hm3$ACTB+1) )
firstfewmonths_hm3$Spoly<- log10( (firstfewmonths_hm3$Spoly+1) / (firstfewmonths_hm3$ACTB+1) )
firstfewmonths_hm3$Srbd<- log10( (firstfewmonths_hm3$Srbd+1) / (firstfewmonths_hm3$ACTB+1) )
firstfewmonths_hm3$ACTB<- log10( firstfewmonths_hm3$ACTB+1 )
#make hm
yorder = c("ACTB",	"Rdrp",	"Spoly",	"Srbd")

firstfewmonths_hm3$sample<-rownames(firstfewmonths_hm3)
firstfewmonths_hm3_long<-melt(firstfewmonths_hm3, id.vars = c("sample"))
colnames(firstfewmonths_hm3_long)<-c("sample", "gene","counts");

firstfewmonths_hm3_long$ordered_y<-firstfewmonths_hm3_long$gene
firstfewmonths_hm3_long$ordered_y <- factor(firstfewmonths_hm3_long$ordered_y, levels = yorder)
firstfewmonths_hm3_long$ordered_x<-firstfewmonths_hm3_long$sample
firstfewmonths_hm3_long$ordered_x <- factor(firstfewmonths_hm3_long$ordered_x, levels = c(orderedbyzscore, as.factor(pilotnegs$sample))) #orderedbydate #x_order_names

##in this vers the samples on the right are the Negatives.
heatmap4<-ggplot(firstfewmonths_hm3_long, aes(x=ordered_x, y=ordered_y, fill = (counts) )) +
  ggtitle("Pilot Set: Log10[(counts+1)/(ACTB+1)], ordered by RBD Zscore") +
  geom_tile() + xlab("Sample") + ylab("Gene") + labs(fill="Zscore Log10(counts+1)") + theme_bw() +
  scale_fill_gradientn(colors=c(low = 'skyblue4', mid = 'wheat2', high = 'orangered3'), limits=c(-5, 5), breaks=c(-5, -2.5, 0, 2.5, 5) )+
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      axis.text.x = element_text(size = 1, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 10)  )
heatmap4

#this is vers in manuscript
pdf("pilotsamples_withnegatives_countHeatmaps_log10actb-norm_byRBDscore_v5_Apr2023.pdf", height = 4, width = 20)
heatmap4
dev.off()
