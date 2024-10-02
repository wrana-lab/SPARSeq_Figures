


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
