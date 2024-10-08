#This is an example using our Delta subset of barcoded data, which shows how we aggregate the results of running our quasispecies 
#pipeline at 7 different numerical read percent cutoffs to ultimately fit a model to the lowest three cutoffs and use that
#to determine an x-intercept value, which we use as a conservative cutoff to run a final pass of the pipeline to detect putative quasispecies (pQS).

#Inputs for this example are provided:

library(reshape2)
library(gridExtra)

##the sparseq quasispecies pipeline needs to be run with the percent cutoff setting changed for each of these outputs. 
#import each output with modified names 
bc1_001<-read.csv("srbd_aligned_list_bc1_0.01.txt")
head(bc1_001)
bc2_001<-read.csv("srbd_aligned_list_bc2_0.01.txt")
head(bc2_001)
bc1_001_list<-read.csv("countoutputlist_bc1_0.01.csv") #this is always outputted from the QS pipeline 
bc2_001_list<-read.csv("countoutputlist_bc2_0.01.csv") #
bc1_001_keeplist<-bc1_001_list[bc1_001_list$total_SRBD_count > 32000,] #32k cutoff as in Sup Fig 12A
bc2_001_keeplist<-bc2_001_list[bc2_001_list$total_SRBD_count > 32000,]
nrow(bc1_001_list); nrow(bc1_001_keeplist) 
nrow(bc2_001_keeplist) 
##pass number is same for all % versions because it's a flat 32k minimum 

bc1_002<-read.csv("srbd_aligned_list_bc1_0.02.txt")
bc2_002<-read.csv("srbd_aligned_list_bc2_0.02.txt")

bc1_005<-read.csv("srbd_aligned_list_bc1_0.05.txt")
bc2_005<-read.csv("srbd_aligned_list_bc2_0.05.txt")

bc1_01<-read.csv("srbd_aligned_list_bc1_0.1.txt")
bc2_01<-read.csv("srbd_aligned_list_bc2_0.1.txt")

bc1_02<-read.csv("srbd_aligned_list_bc1_0.2.txt")
bc2_02<-read.csv("srbd_aligned_list_bc2_0.2.txt")

bc1_05<-read.csv("srbd_aligned_list_bc1_0.5.txt")
bc2_05<-read.csv("srbd_aligned_list_bc2_0.5.txt")

bc1_1<-read.csv("srbd_aligned_list_bc1_1.txt") 
bc2_1<-read.csv("srbd_aligned_list_bc2_1.txt")

####graph of counts cutoffs#####
bc1_001_list<-bc1_001_list[order(bc1_001_list$total_SRBD_count, decreasing = T),]
bc2_001_list<-bc2_001_list[order(bc2_001_list$total_SRBD_count, decreasing = T),]
bc1_001_list$sample_F<-factor(bc1_001_list$sample, levels=bc1_001_list$sample)
bc2_001_list$sample_F<-factor(bc2_001_list$sample, levels=bc2_001_list$sample)

sampleplot1a<-ggplot(bc1_001_list, aes(x = sample_F, y=total_SRBD_count)) + scale_y_log10(limits = c(1, 500000)) +
  geom_point() + ggtitle("Srbm Counts BC1") + geom_hline(yintercept=32000) +
  ylab("SRBD count") + xlab("Sample") + theme_bw() + 
  theme(axis.title = element_text(size=18), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5, size = 20),
      panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 14))
sampleplot1a

sampleplot2a<-ggplot(bc2_001_list, aes(x = sample_F, y=total_SRBD_count)) +  scale_y_log10(limits = c(1, 500000)) +
  geom_point() + ggtitle("Srbm Counts BC2") + geom_hline(yintercept=32000) + 
  ylab("SRBD count") + xlab("Sample") + theme_bw() + 
  theme(axis.title = element_text(size=18), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5, size = 20),
      panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 14))
sampleplot2a

pdf("example_Srbm_counts_Deltadata.pdf", width = 12, height = 5)
grid.arrange(sampleplot1a, sampleplot2a, ncol=2)
dev.off()

####aggregate and fit linear model to the lowest three cutoffs####
#drop refseq from top row as we only used it for aligning and visualizing the intermediate outputs 
bc1_001<-bc1_001[2:nrow(bc1_001),]
bc2_001<-bc2_001[2:nrow(bc2_001),]
bc1_002<-bc1_002[2:nrow(bc1_002),]
bc2_002<-bc2_002[2:nrow(bc2_002),]
bc1_005<-bc1_005[2:nrow(bc1_005),]
bc2_005<-bc2_005[2:nrow(bc2_005),]
bc1_01<-bc1_01[2:nrow(bc1_01),]
bc2_01<-bc2_01[2:nrow(bc2_01),]
bc1_02<-bc1_02[2:nrow(bc1_02),]
bc2_02<-bc2_02[2:nrow(bc2_02),]
bc1_05<-bc1_05[2:nrow(bc1_05),]
bc2_05<-bc2_05[2:nrow(bc2_05),]
bc1_1<-bc1_1[2:nrow(bc1_1),]
bc2_1<-bc2_1[2:nrow(bc2_1),]

#organize tables, adjust sample ID and drop gaps from sequences 
bc1_001$sampleID<-substr(bc1_001$sample, 1, 9)
bc1_001$seq_chars<-gsub("-", "", bc1_001$sequence); bc1_001$seq_chars<-gsub("C$", "", bc1_001$seq_chars)
bc2_001$sampleID<-substr(bc2_001$sample, 1, 9)
bc2_001$seq_chars<-gsub("-", "", bc2_001$sequence); bc2_001$seq_chars<-gsub("C$", "", bc2_001$seq_chars)
bc1_001_sm<-bc1_001[,3:4]; bc2_001_sm<-bc2_001[,3:4]

bc1_002$sampleID<-substr(bc1_002$sample, 1, 9)
bc1_002$seq_chars<-gsub("-", "", bc1_002$sequence); bc1_002$seq_chars<-gsub("C$", "", bc1_002$seq_chars)
bc2_002$sampleID<-substr(bc2_002$sample, 1, 9)
bc2_002$seq_chars<-gsub("-", "", bc2_002$sequence); bc2_002$seq_chars<-gsub("C$", "", bc2_002$seq_chars)
bc1_002_sm<-bc1_002[,3:4]; bc2_002_sm<-bc2_002[,3:4]

bc1_005$sampleID<-substr(bc1_005$sample, 1, 9)
bc1_005$seq_chars<-gsub("-", "", bc1_005$sequence); bc1_005$seq_chars<-gsub("C$", "", bc1_005$seq_chars)
bc2_005$sampleID<-substr(bc2_005$sample, 1, 9)
bc2_005$seq_chars<-gsub("-", "", bc2_005$sequence); bc2_005$seq_chars<-gsub("C$", "", bc2_005$seq_chars)
bc1_005_sm<-bc1_005[,3:4]; bc2_005_sm<-bc2_005[,3:4]

bc1_01$sampleID<-substr(bc1_01$sample, 1, 9)
bc1_01$seq_chars<-gsub("-", "", bc1_01$sequence); bc1_01$seq_chars<-gsub("C$", "", bc1_01$seq_chars)
bc2_01$sampleID<-substr(bc2_01$sample, 1, 9)
bc2_01$seq_chars<-gsub("-", "", bc2_01$sequence); bc2_01$seq_chars<-gsub("C$", "", bc2_01$seq_chars)
bc1_01_sm<-bc1_01[,3:4]; bc2_01_sm<-bc2_01[,3:4]

bc1_02$sampleID<-substr(bc1_02$sample, 1, 9)
bc1_02$seq_chars<-gsub("-", "", bc1_02$sequence); bc1_02$seq_chars<-gsub("C$", "", bc1_02$seq_chars)
bc2_02$sampleID<-substr(bc2_02$sample, 1, 9)
bc2_02$seq_chars<-gsub("-", "", bc2_02$sequence); bc2_02$seq_chars<-gsub("C$", "", bc2_02$seq_chars)
bc1_02_sm<-bc1_02[,3:4]; bc2_02_sm<-bc2_02[,3:4]

bc1_05$sampleID<-substr(bc1_05$sample, 1, 9)
bc1_05$seq_chars<-gsub("-", "", bc1_05$sequence); bc1_05$seq_chars<-gsub("C$", "", bc1_05$seq_chars)
bc2_05$sampleID<-substr(bc2_05$sample, 1, 9)
bc2_05$seq_chars<-gsub("-", "", bc2_05$sequence); bc2_05$seq_chars<-gsub("C$", "", bc2_05$seq_chars)
bc1_05_sm<-bc1_05[,3:4]; bc2_05_sm<-bc2_05[,3:4]

bc1_1$sampleID<-substr(bc1_1$sample, 1, 9)
bc1_1$seq_chars<-gsub("-", "", bc1_1$sequence); bc1_1$seq_chars<-gsub("C$", "", bc1_1$seq_chars)
bc2_1$sampleID<-substr(bc2_1$sample, 1, 9)
bc2_1$seq_chars<-gsub("-", "", bc2_1$sequence); bc2_1$seq_chars<-gsub("C$", "", bc2_1$seq_chars)
bc1_1_sm<-bc1_1[,3:4]; bc2_1_sm<-bc2_1[,3:4]

##need to add plate info/percent cutoffs for later graph
bc1_001_sm$percent<-0.01
bc2_001_sm$percent<-0.01
bc1_002_sm$percent<-0.02
bc2_002_sm$percent<-0.02
bc1_005_sm$percent<-0.05
bc2_005_sm$percent<-0.05
bc1_01_sm$percent<-0.1
bc2_01_sm$percent<-0.1
bc1_02_sm$percent<-0.2
bc2_02_sm$percent<-0.2
bc1_05_sm$percent<-0.5
bc2_05_sm$percent<-0.5
bc1_1_sm$percent<-1
bc2_1_sm$percent<-1

bc1s<-rbind(bc1_001_sm, bc1_002_sm, bc1_005_sm, bc1_01_sm, bc1_02_sm, bc1_05_sm, bc1_1_sm)
bc2s<-rbind(bc2_001_sm, bc2_002_sm, bc2_005_sm, bc2_01_sm, bc2_02_sm, bc2_05_sm, bc2_1_sm)


##set up graphs showing the number of unique sequences detected in bc1/bc2 at each percent cutoff 
bc1_count_agg<-as.data.frame(bc1s %>% group_by(percent) %>%  summarise(qscountsbc1 = n_distinct(seq_chars)) )
bc2_count_agg<-as.data.frame(bc2s %>% group_by(percent) %>%  summarise(qscountsbc2 = n_distinct(seq_chars)) )

bc1_count_agg<-as.data.frame(bc1s %>% group_by(percent) %>%  summarise(qscountsbc1 = n()) )
bc2_count_agg<-as.data.frame(bc2s %>% group_by(percent) %>%  summarise(qscountsbc2 = n()) )

run_summary_both<-merge(bc1_count_agg,bc2_count_agg,by="percent")
run_summary_both<-melt(run_summary_both, id.vars=c("percent"))
#set exp to just be BC1 vs BC2 as we don't colour by the original clinical run here anymore 
run_summary_both$Exp<-"BC1"
run_summary_both[run_summary_both$variable == "qscountsbc2",]$Exp<-"BC2"


qsplot_facet<-ggplot(run_summary_both, aes(x = percent, y=value, colour = Exp)) + facet_wrap(vars(Exp)) +
  geom_point(size=2,alpha=1) + ggtitle("Delta QS Summary:\n32k min. SRbm counts") + 
  scale_y_continuous(breaks = c(seq(0,30000,5000)), limits=c(0,30000)) + ylab("QS Number") + xlab("Percent Cutoff") + theme_bw() +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(hjust = 0.5, size = 10), axis.text.y = element_text(hjust = 1, size = 10))
qsplot_facet

pdf("delta_qs_percent_cutoffs_normal_yaxis.pdf", width = 8, height = 5)
qsplot_facet
dev.off()


qsplot_facetB<-ggplot(run_summary_both, aes(x = percent, y=log2(value+1), colour = Exp)) + facet_wrap(vars(Exp)) +  geom_point(size=2,alpha=1) + 
  ggtitle("Delta QS Summary:\n32k min. SRBD counts") + ylab("QS Number as Log2(#+1)") + xlab("Percent Cutoff") + theme_bw() +xlim(0,1) +  ylim(0,20) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(hjust = 0.5, size = 10), axis.text.y = element_text(hjust = 1, size = 10))
qsplot_facetB

pdf("delta_qs_percent_cutoffs_log_yaxis.pdf", width = 10, height = 5)
qsplot_facetB
dev.off()


###nonsplit version
run_summary_plate_both<-merge(bc1_count_agg,bc2_count_agg, by="percent")
run_summary_plate_both<-melt(run_summary_plate_both, id.vars=c("percent"))
colnames(run_summary_plate_both)<-c( "percent", "variable", "value")
run_summary_plate_both$BC<-""
run_summary_plate_both[run_summary_plate_both$variable == "qscountsbc1",]$BC<-"BC1"
run_summary_plate_both[run_summary_plate_both$variable == "qscountsbc2",]$BC<-"BC2"

###get x-int of model of first 3 cutoffs- we found this provides a conservative cutoff that gives a strict results list with little background
bc1_delta_bc1_small<-subset(run_summary_plate_both,run_summary_plate_both$BC == "BC1" & run_summary_plate_both$percent %in% c(0.01,0.02,0.05))
bc1_delta_bc2_small<-subset(run_summary_plate_both,run_summary_plate_both$BC == "BC2" & run_summary_plate_both$percent %in% c(0.01,0.02,0.05))
bc1_delta_bc1_small$log2counts<-log2(bc1_delta_bc1_small$value)
bc1_delta_bc2_small$log2counts<-log2(bc1_delta_bc2_small$value)

###apply linear model - quickly visualize example
ggplot(bc1_delta_bc1_small, aes(x = percent, y=(log2counts))) + geom_point(size=2,alpha=0.5) + ggtitle("Delta QS Summary: 32k min. SRbm counts") + 
   scale_y_continuous(limits=c(0,20)) + xlim(0,0.2) + ylab("log2(QS Number)") + xlab("Percent Cutoff") + theme_bw() + 
   theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
       panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
       panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
       axis.text.x = element_text(hjust = 0.5, size = 10), axis.text.y = element_text(hjust = 1, size = 10)) +
   stat_smooth(method = "lm", formula = y ~ x, se=F)


estimate_bc1_delta_bc1_small <- lm((log2counts) ~ percent, data=bc1_delta_bc1_small)
summary(estimate_bc1_delta_bc1_small)
# Residuals:
#        1        2        3 
#  0.18583 -0.24778  0.06194 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  14.9553     0.3393  44.079   0.0144 *
# percent     -59.5519    10.7291  -5.551   0.1135  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3159 on 1 degrees of freedom
# Multiple R-squared:  0.9686,	Adjusted R-squared:  0.9371 
# F-statistic: 30.81 on 1 and 1 DF,  p-value: 0.1135

#y = mx + b
# y int is 14.9553
#slope is -59.5519 
# 0 = -59.5519 * x + 14.3991
# -14.9553 /  -59.5519 = 0.2511305
##so x int for cutoff is  0.2511


bc1_delta_bc2_small <- lm((log2counts) ~ percent, data=bc1_delta_bc2_small)
summary(bc1_delta_bc2_small)
# Residuals:
#        8        9       10 
#  0.20718 -0.27624  0.06906 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  15.0271     0.3783  39.727    0.016 *
# percent     -74.1307    11.9616  -6.197    0.102  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3521 on 1 degrees of freedom
# Multiple R-squared:  0.9746,	Adjusted R-squared:  0.9492 
# F-statistic: 38.41 on 1 and 1 DF,  p-value: 0.1018

#y = mx + b
# y int is 15.0271
#slope is -74.1307 
# 0 = -74.1307 * x + 15.0271
# -15.0271 /  -74.1307 = 0.2027109
##so x int for cutoff = 0.2027

#generate plot showing the x-int cutoffs 
qsplot_facetB2<-ggplot(run_summary_both, aes(x = percent, y=log2(value+1), colour = Exp)) + facet_wrap(vars(Exp)) + 
  geom_abline(data = subset(run_summary_both, Exp=="BC1"), aes(slope = -59.5519 , intercept = 14.9553))+ 
  geom_abline(data = subset(run_summary_both, Exp=="BC2"), aes(slope = -74.1307, intercept = 15.0271))+ 
   geom_point(size=2,alpha=1) + ggtitle("139-141 Delta SRBM QS Summary: 10k min. SRBM counts") + 
    ylab("QS Number as Log2(#+1)") + xlab("Percent Cutoff") + theme_bw() +  xlim(0,1) + ylim(0, 20) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5), 
      panel.border = element_rect(colour = "black", fill=NA, size=1),  plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(hjust = 0.5, size = 12), axis.text.y = element_text(hjust = 1, size = 12))
qsplot_facetB2
pdf("SRBM_Delta_QS_modelfitwithline.pdf", width = 7,height =5)
qsplot_facetB2
dev.off()
