###supp 3 and 4 use run 28-31 data so we do them together####

#We compare the counts tables for run 28, 29, 30, and 31 since we assessed them with our V1 primers and V2 primers
###to show that the change to V2 primers is ok

library(openxlsx)
library(ggplot2)
options(scipen=999)


v1cl28<-read.xlsx("runcl28_counts.xlsx"); v1cl28$run<-"28"
v1cl29<-read.xlsx("runcl29_v1_countTable.xlsx"); v1cl29$run<-"29"
v1cl30<-read.xlsx("runcl30_countTable.xlsx"); v1cl30$run<-"30"
v1cl31<-read.xlsx("runcl31_countTable.xlsx"); v1cl31$run<-"31"

v2cl28<-read.xlsx("runcl28_v2_CountTable.xlsx"); v2cl28$run<-"28"
v2cl29<-read.xlsx("runcl29_v2_CountTable.xlsx"); v2cl29$run<-"29"
v2cl30<-read.xlsx("runcl30_v2_CountTable.xlsx"); v2cl30$run<-"30"
v2cl31<-read.xlsx("runcl31_v2_CountTable.xlsx"); v2cl31$run<-"31"

tokeep1<-c("sample","ACTB","Rdrp","Spoly","Srbd","total.raw.reads","run")
v1cl28<-v1cl28[,tokeep1]; v1cl29<-v1cl29[,tokeep1]
v1cl30<-v1cl30[,tokeep1]; v1cl31<-v1cl31[,tokeep1]
allv1<-rbind(v1cl28, v1cl29, v1cl30, v1cl31)
allv1$Srbd_v2<-0
allv1<-allv1[,c("sample","ACTB","Rdrp","Spoly","Srbd","Srbd_v2","total.raw.reads","run")]
colnames(allv1)<-c("sample","v1.ACTB","v1.Rdrp","v1.Spoly","v1.Srbd", "v1.Srbd_v2","v1.total.raw.reads","v1.run")

tokeep2<-c("sample","ACTB","Rdrp","Spoly","Srbd","Srbd_v2","total.raw.reads","run")
v2cl28<-v2cl28[,tokeep2]; v2cl29<-v2cl29[,tokeep2]
v2cl30<-v2cl30[,tokeep2]; v2cl31<-v2cl31[,tokeep2]
allv2<-rbind(v2cl28, v2cl29, v2cl30, v2cl31)
colnames(allv2)<-c("sample","v2.ACTB","v2.Rdrp","v2.Spoly","v2.Srbd", "v2.Srbd_v2","v2.total.raw.reads","v2.run")

#make sample_run for merging to work
allv1$sample_run<-paste(allv1$sample, allv1$v1.run, sep = "_")
allv2$sample_run<-paste(allv2$sample, allv2$v2.run, sep = "_")

allv1v2<-merge(allv1,allv2, by = "sample_run", all = T)
#anything with an X is a sample, else is control
allv1v2$type<-"Control"
allv1v2[grep("X", allv1v2$sample_run, invert = F),]$type<-"Sample"
nrow(allv1v2) #146 

##drop any that aren't either a control or found in actual dataset
allv1v2<-subset(allv1v2, (allv1v2$type == "Control") | (allv1v2$sample.x %in% re_filtered_maindataset$sample))
nrow(allv1v2) #122

#create scatterplot of srbd vs srbd_v2, colour by run and shape by sample or control

#all data
cor.test(allv1v2$v1.Srbd, allv1v2$v2.Srbd_v2, method = 'spearman')
#rho = 0.87, pval < 0.00000000000000022
# 0.8671036^2 R2 is 0.75

#control only
allv1v2.controls<-subset(allv1v2, allv1v2$type == "Control")
cor.test(allv1v2.controls$v1.Srbd, allv1v2.controls$v2.Srbd_v2, method = 'spearman')
#rho = 0.427481, pval = 0.02327

#sample only > this is the kept vers 
allv1v2.samples<-subset(allv1v2, allv1v2$type == "Sample")
cor.test(allv1v2.samples$v1.Srbd, allv1v2.samples$v2.Srbd_v2, method = 'spearman')
#rho = 0.8161926, p-value < 0.00000000000000022
#R2 = 0.8161926^2 == 0.666

#allv1v2<-subset(allv1v2, allv1v2$type != "Control")
nrow(allv1v2.samples) #94

plotSrbd.Srbdv2<-ggplot(allv1v2.samples, aes(x=v1.Srbd, y=v2.Srbd_v2, colour = v1.run)) +
  geom_point(alpha = 0.85) + ggtitle("Srbd Counts, Clinical Runs 28-31\n[Spearman cor. = 0.82 for samples only,\nR2 = 0.67]") +
  xlab("V1 Srbd") + ylab("V2 Srbd") + theme_bw() + scale_x_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000), limits = c(1,1000000)) + 
  scale_y_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000,10000000), limits = c(1, 10000000)) + coord_equal() +  theme_bw() +  labs(colour="Run ID") + 
  theme(axis.title = element_text(size=16), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5),  panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
plotSrbd.Srbdv2


pdf("Srbd-Srbdv2_Comparison_SamplesOnly_v25data.pdf", width=6, height = 6)
plotSrbd.Srbdv2
dev.off()

###repeat above final version but for RDRP, SPBS, ACTB

#sample only
allv1v2.samples<-subset(allv1v2, allv1v2$type == "Sample")
cor.test(allv1v2.samples$v1.Rdrp, allv1v2.samples$v2.Rdrp, method = 'spearman')
#r2 = 0.8559228^2 == 0.73

table(allv1v2$type)

plotRdrp<-ggplot(allv1v2.samples, aes(x=v1.Rdrp, y=v2.Rdrp, colour = v1.run)) +
  geom_point(alpha = 0.85) + ggtitle("Rdrp Counts, Clinical Runs 28-31\n[Spearman cor. = 0.86 for samples only,\nR2 = 0.73]") +
  xlab("V1 Rdrp") + ylab("V2 Rdrp") + theme_bw() + scale_x_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000, 10000000), limits = c(1,10000000)) + 
  scale_y_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000,10000000), limits = c(1, 10000000)) + coord_equal() +  theme_bw() +  labs(colour="Run ID") +
  theme(axis.title = element_text(size=16), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5),      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
plotRdrp
pdf("RDRP_Comparison_SamplesOnly_v25data.pdf", width=6, height = 6)
plotRdrp
dev.off()

cor.test(allv1v2.samples$v1.Spoly, allv1v2.samples$v2.Spoly, method = 'spearman')
#0.810772 
#r2 == 0.810772^2 = 0.657

plotSPBS<-ggplot(allv1v2.samples, aes(x=v1.Spoly, y=v2.Spoly, colour = v1.run)) +
  geom_point(alpha = 0.85) + ggtitle("Spbs Counts, Clinical Runs 28-31\n[Spearman cor. = 0.81 for samples only,\nR2 = 0.657]") +
  xlab("V1 Spbs") + ylab("V2 Spbs") + theme_bw() + scale_x_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000, 10000000), limits = c(1,10000000)) + 
  scale_y_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000,10000000), limits = c(1, 10000000)) + coord_equal() +  theme_bw() +  labs(colour="Run ID") + 
  theme(axis.title = element_text(size=16), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
plotSPBS

pdf("SPBS_Comparison_SamplesOnly_v25data.pdf", width=6, height = 6)
plotSPBS
dev.off()


cor.test(allv1v2.samples$v1.ACTB, allv1v2.samples$v2.ACTB, method = 'spearman')
# r2 = 0.8914095^2 == 0.795

plotACTB<-ggplot(allv1v2.samples, aes(x=v1.ACTB, y=v2.ACTB, colour = v1.run)) +
  geom_point(alpha = 0.85) + ggtitle("ACTB Counts, Clinical Runs 28-31\n[Spearman cor. = 0.88 for samples only,\nR2 = 0.795]") +
  xlab("V1 ACTB") + ylab("V2 ACTB") + theme_bw() + scale_x_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000, 10000000), limits = c(1,10000000)) + 
  scale_y_log10(breaks = c( 1,10,100,1000,10000,100000, 1000000,10000000), limits = c(1, 10000000)) + coord_equal() +  theme_bw() +
  labs(colour="Run ID") +
  theme(axis.title = element_text(size=16), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1, size = 10))
plotACTB

pdf("ACTB_Comparison_SamplesOnly_v25.pdf", width=6, height = 6)
plotACTB
dev.off()
