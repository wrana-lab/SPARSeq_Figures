#Script to aggregate results of individual pQS pipeline runs using customized cutoffs as described in previous scripts.

#load all outputs
finalwt27seqcounts<-read.xlsx("srbd_outputs/wt27_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized.xlsx")
finalwt31seqcounts<-read.xlsx("srbd_outputs/wt31_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized.xlsx")
finalalpha31seqcounts<-read.xlsx("srbd_outputs/alpha31_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized.xlsx")
finaldeltaseqcounts<-read.xlsx("srbd_outputs/delta_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized.xlsx")
finalomicronseqcounts<-read.xlsx("srbd_outputs/omicron_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized.xlsx")

finalwt27samplesseq<-read.xlsx("srbd_outputs/wt27_pQS_sample_sequence_list.xlsx")
finalwt31samplesseq<-read.xlsx("srbd_outputs/wt31_pQS_sample_sequence_list.xlsx")
finalalpha31samplesseq<-read.xlsx("srbd_outputs/alpha31_pQS_sample_sequence_list.xlsx")
finaldeltasamplesseq<-read.xlsx("srbd_outputs/delta_pQS_sample_sequence_list.xlsx")
finalomicronsamplesseq<-read.xlsx("srbd_outputs/omicron_pQS_sample_sequence_list.xlsx")

finalwt27seqcounts$variant<-"WT"
finalwt31seqcounts$variant<-"WT"
finalalpha31seqcounts$variant<-"Alpha"
finaldeltaseqcounts$variant<-"Delta"
finalomicronseqcounts$variant<-"Omicron"

finalwt27samplesseq$variant<-"WT"; finalwt27samplesseq$run<-"27"
finalwt31samplesseq$variant<-"WT"; finalwt31samplesseq$run<-"31"
finalalpha31samplesseq$variant<-"Alpha"; finalalpha31samplesseq$run<-"31"
finaldeltasamplesseq$variant<-"Delta"; finaldeltasamplesseq$run<-"139-141"
finalomicronsamplesseq$variant<-"Omicron"; finalomicronsamplesseq$run<-"139-141"

srbdoverall<-rbind(finalwt27samplesseq, finalwt31samplesseq, finalalpha31samplesseq, finaldeltasamplesseq, finalomicronsamplesseq)


#We found that the sample producing the single pQS for SPBS for run 27 WT would have been annotated WT+ 
#based on results from barcode runs so decided to exclude it as an artefact
#therefore run 27 SPBS is not included for downstream analysis 

#finalwt27seqcounts_spbs<-read.xlsx("spbs_outputs/wt27_pQS_cutoffs_common_sampleIDsandsequences_finalized_spbs.xlsx")
finalwt31seqcounts_spbs<-read.xlsx("spbs_outputs/wt31_pQS_cutoffs_common_sampleIDsandsequences_finalized_spbs.xlsx")
finalalpha31seqcounts_spbs<-read.xlsx("spbs_outputs/alpha31_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized_spbs.xlsx")
finaldeltaseqcounts_spbs<-read.xlsx("spbs_outputs/delta_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized_spbs.xlsx")
finalomicronseqcounts_spbs<-read.xlsx("spbs_outputs/omicron_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized_spbs.xlsx")

#finalwt27samplesseq_spbs<-read.xlsx("spbs_outputs/wt27_pQS_sample_sequence_list_spbs.xlsx")
finalwt31samplesseq_spbs<-read.xlsx("spbs_outputs/wt31_pQS_sample_sequence_list_spbs.xlsx")
finalalpha31samplesseq_spbs<-read.xlsx("spbs_outputs/alpha31_pQS_sample_sequence_list_spbs.xlsx")
finaldeltasamplesseq_spbs<-read.xlsx("spbs_outputs/delta_pQS_sample_sequence_list.xlsx")
finalomicronsamplesseq_spbs<-read.xlsx("spbs_outputs/omicron_pQS_sample_sequence_list_spbs.xlsx")
#for omicron we had adjusted spbs because the beginning of the amplicon had a mutation or something so this file was _adjusted until i renamed it here 

#finalwt27seqcounts_spbs$variant<-"wt"
finalwt31seqcounts_spbs$variant<-"WT"
finalalpha31seqcounts_spbs$variant<-"Alpha"
finaldeltaseqcounts_spbs$variant<-"Delta"
finalomicronseqcounts_spbs$variant<-"Omicron"

#finalwt27samplesseq_spbs$variant<-"WT"; finalwt27samplesseq_spbs$run<-"27"
finalwt31samplesseq_spbs$variant<-"WT"; finalwt31samplesseq_spbs$run<-"31"
finalalpha31samplesseq_spbs$variant<-"Alpha"; finalalpha31samplesseq_spbs$run<-"31"
finaldeltasamplesseq_spbs$variant<-"Delta"; finaldeltasamplesseq_spbs$run<-"139-141"
finalomicronsamplesseq_spbs$variant<-"Omicron"; finalomicronsamplesseq_spbs$run<-"139-141"

#finalwt27samplesseq_spbs,  is excluded, see above
spbsoverall<-rbind(finalwt31samplesseq_spbs, finalalpha31samplesseq_spbs, finaldeltasamplesseq_spbs, finalomicronsamplesseq_spbs)


head(srbdoverall)
head(spbsoverall)
srbdoverall$sample_seq<-paste(srbdoverall$sampleID, srbdoverall$seq_chars, sep = "_")
spbsoverall$sample_seq<-paste(spbsoverall$sampleID, spbsoverall$seq_chars, sep = "_")

##read in all the seq percent lists...
srbd_wt27bc1percents<-read.csv("srbd_outputs/seq_percents_srbd_run27WT_bc1.csv", header = F)
srbd_wt27bc2percents<-read.csv("srbd_outputs/seq_percents_srbd_run27WT_bc2.csv", header = F)
srbd_wt27bc1percents<-srbd_wt27bc1percents[,c(3,4,6)]
srbd_wt27bc2percents<-srbd_wt27bc2percents[,c(3,4,6)]

srbd_wt31bc1percents<-read.csv("srbd_outputs/seq_percents_srbd_run31WT_bc1.csv", header = F)
srbd_wt31bc2percents<-read.csv("srbd_outputs/seq_percents_srbd_run31WT_bc2.csv", header = F)
srbd_wt31bc1percents<-srbd_wt31bc1percents[,c(3,4,6)]
srbd_wt31bc2percents<-srbd_wt31bc2percents[,c(3,4,6)]

srbd_alpha31bc1percents<-read.csv("srbd_outputs/seq_percents_srbd_run31Alpha_bc1.csv", header = F)
srbd_alpha31bc2percents<-read.csv("srbd_outputs/seq_percents_srbd_run31Alpha_bc2.csv", header = F)
srbd_alpha31bc1percents<-srbd_alpha31bc1percents[,c(3,4,6)]
srbd_alpha31bc2percents<-srbd_alpha31bc2percents[,c(3,4,6)]

srbd_deltabc1percents<-read.csv("srbd_outputs/seq_percents_srbd_runDelta_bc1.csv", header = F)
srbd_deltabc2percents<-read.csv("srbd_outputs/seq_percents_srbd_runDelta_bc2.csv", header = F)
srbd_deltabc1percents<-srbd_deltabc1percents[,c(3,4,6)]
srbd_deltabc2percents<-srbd_deltabc2percents[,c(3,4,6)]

srbd_omicronbc1percents<-read.csv("srbd_outputs/seq_percents_srbd_runOmicron_bc1.csv", header = F)
srbd_omicronbc2percents<-read.csv("srbd_outputs/seq_percents_srbd_runOmicron_bc2.csv", header = F)
srbd_omicronbc1percents<-srbd_omicronbc1percents[,c(3,4,6)]
srbd_omicronbc2percents<-srbd_omicronbc2percents[,c(3,4,6)]

#spbs
spbs_wt27bc1percents<-read.csv("spbs_outputs/seq_percents_spbs_run27WT_bc1.csv", header = F)
spbs_wt27bc2percents<-read.csv("spbs_outputs/seq_percents_spbs_run27WT_bc2.csv", header = F)
spbs_wt27bc1percents<-spbs_wt27bc1percents[,c(3,4,6)]
spbs_wt27bc2percents<-spbs_wt27bc2percents[,c(3,4,6)]

spbs_wt31bc1percents<-read.csv("spbs_outputs/seq_percents_spbs_run31WT_bc1.csv", header = F)
spbs_wt31bc2percents<-read.csv("spbs_outputs/seq_percents_spbs_run31WT_bc2.csv", header = F)
spbs_wt31bc1percents<-spbs_wt31bc1percents[,c(3,4,6)]
spbs_wt31bc2percents<-spbs_wt31bc2percents[,c(3,4,6)]

spbs_alpha31bc1percents<-read.csv("spbs_outputs/seq_percents_spbs_run31Alpha_bc1.csv", header = F)
spbs_alpha31bc2percents<-read.csv("spbs_outputs/seq_percents_spbs_run31Alpha_bc2.csv", header = F)
spbs_alpha31bc1percents<-spbs_alpha31bc1percents[,c(3,4,6)]
spbs_alpha31bc2percents<-spbs_alpha31bc2percents[,c(3,4,6)]

spbs_deltabc1percents<-read.csv("spbs_outputs/seq_percents_spbs_Delta_bc1.csv", header = F)
spbs_deltabc2percents<-read.csv("spbs_outputs/seq_percents_spbs_Delta_bc2.csv", header = F)
spbs_deltabc1percents<-spbs_deltabc1percents[,c(3,4,6)]
spbs_deltabc2percents<-spbs_deltabc2percents[,c(3,4,6)]

spbs_omicronbc1percents<-read.csv("spbs_outputs/seq_percents_spbs_Omicron_bc1.csv", header = F)
spbs_omicronbc2percents<-read.csv("spbs_outputs/seq_percents_spbs_Omicron_bc2.csv", header = F)
spbs_omicronbc1percents<-spbs_omicronbc1percents[,c(3,4,6)]
spbs_omicronbc2percents<-spbs_omicronbc2percents[,c(3,4,6)]

##add sampleID-seq column to all of these
srbd_wt27bc1percents$sample_seq<-paste(srbd_wt27bc1percents$V6, srbd_wt27bc1percents$V3, sep = "_")
srbd_wt27bc2percents$sample_seq<-paste(srbd_wt27bc2percents$V6, srbd_wt27bc2percents$V3, sep = "_")

srbd_wt31bc1percents$sample_seq<-paste(srbd_wt31bc1percents$V6, srbd_wt31bc1percents$V3, sep = "_")
srbd_wt31bc2percents$sample_seq<-paste(srbd_wt31bc2percents$V6, srbd_wt31bc2percents$V3, sep = "_")

srbd_alpha31bc1percents$sample_seq<-paste(srbd_alpha31bc1percents$V6, srbd_alpha31bc1percents$V3, sep = "_")
srbd_alpha31bc2percents$sample_seq<-paste(srbd_alpha31bc2percents$V6, srbd_alpha31bc2percents$V3, sep = "_")

srbd_deltabc1percents$sample_seq<-paste(srbd_deltabc1percents$V6, srbd_deltabc1percents$V3, sep = "_")
srbd_deltabc2percents$sample_seq<-paste(srbd_deltabc2percents$V6, srbd_deltabc2percents$V3, sep = "_")

srbd_omicronbc1percents$sample_seq<-paste(srbd_omicronbc1percents$V6, srbd_omicronbc1percents$V3, sep = "_")
srbd_omicronbc2percents$sample_seq<-paste(srbd_omicronbc2percents$V6, srbd_omicronbc2percents$V3, sep = "_")


spbs_wt27bc1percents$sample_seq<-paste(spbs_wt27bc1percents$V6, spbs_wt27bc1percents$V3, sep = "_")
spbs_wt27bc2percents$sample_seq<-paste(spbs_wt27bc2percents$V6, spbs_wt27bc2percents$V3, sep = "_")

spbs_wt31bc1percents$sample_seq<-paste(spbs_wt31bc1percents$V6, spbs_wt31bc1percents$V3, sep = "_")
spbs_wt31bc2percents$sample_seq<-paste(spbs_wt31bc2percents$V6, spbs_wt31bc2percents$V3, sep = "_")

spbs_alpha31bc1percents$sample_seq<-paste(spbs_alpha31bc1percents$V6, spbs_alpha31bc1percents$V3, sep = "_")
spbs_alpha31bc2percents$sample_seq<-paste(spbs_alpha31bc2percents$V6, spbs_alpha31bc2percents$V3, sep = "_")

spbs_deltabc1percents$sample_seq<-paste(spbs_deltabc1percents$V6, spbs_deltabc1percents$V3, sep = "_")
spbs_deltabc2percents$sample_seq<-paste(spbs_deltabc2percents$V6, spbs_deltabc2percents$V3, sep = "_")

spbs_omicronbc1percents$sample_seq<-paste(spbs_omicronbc1percents$V6, spbs_omicronbc1percents$V3, sep = "_")
spbs_omicronbc2percents$sample_seq<-paste(spbs_omicronbc2percents$V6, spbs_omicronbc2percents$V3, sep = "_")

#merge each to main srbd or spbs list
srbd_wt27bc1percents<-merge(subset(srbdoverall, srbdoverall$run == "27" & srbdoverall$variant == "WT"), srbd_wt27bc1percents, by = "sample_seq", all= F)
nrow(srbd_wt27bc1percents) #205
table(srbd_wt27bc1percents$V3 == srbd_wt27bc1percents$seq_chars)
srbd_wt27bc1percents<-unique(srbd_wt27bc1percents)

srbd_wt27bc2percents<-merge(subset(srbdoverall, srbdoverall$run == "27" & srbdoverall$variant == "WT"), srbd_wt27bc2percents, by = "sample_seq", all= F)
nrow(srbd_wt27bc2percents) #260
table(srbd_wt27bc2percents$V3 == srbd_wt27bc2percents$seq_chars)
srbd_wt27bc2percents<-unique(srbd_wt27bc2percents)

srbd_wt_27_info<-merge(srbd_wt27bc1percents, srbd_wt27bc2percents, by = "sample_seq", all = F)
head(srbd_wt_27_info,3)
table(srbd_wt_27_info$sampleID.x == srbd_wt_27_info$sampleID.y)
table(srbd_wt_27_info$seq_chars.x == srbd_wt_27_info$seq_chars.y)
srbd_wt_27_info<-srbd_wt_27_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(srbd_wt_27_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")

###we found that W60603525's pQS sequence is just WT which we believe is just an artefact so we removed it	
#ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT
srbd_wt_27_info[srbd_wt_27_info$sample_seq == "W60603525_ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT",]
#remove it 
srbd_wt_27_info<-subset(srbd_wt_27_info, srbd_wt_27_info$sample_seq != "W60603525_ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT")


###srbd 31
srbd_wt31bc1percents<-merge(subset(srbdoverall, srbdoverall$run == "31" & srbdoverall$variant == "WT"), srbd_wt31bc1percents, by = "sample_seq", all= F)
nrow(srbd_wt31bc1percents) #279
table(srbd_wt31bc1percents$V3 == srbd_wt31bc1percents$seq_chars)
srbd_wt31bc1percents<-unique(srbd_wt31bc1percents)

srbd_wt31bc2percents<-merge(subset(srbdoverall, srbdoverall$run == "31" & srbdoverall$variant == "WT"), srbd_wt31bc2percents, by = "sample_seq", all= F)
nrow(srbd_wt31bc2percents) #394now
table(srbd_wt31bc2percents$V3 == srbd_wt31bc2percents$seq_chars)
srbd_wt31bc2percents<-unique(srbd_wt31bc2percents)

srbd_wt_31_info<-merge(srbd_wt31bc1percents, srbd_wt31bc2percents, by = "sample_seq", all = F)
head(srbd_wt_31_info,3)
table(srbd_wt_31_info$sampleID.x == srbd_wt_31_info$sampleID.y)
table(srbd_wt_31_info$seq_chars.x == srbd_wt_31_info$seq_chars.y)
srbd_wt_31_info<-srbd_wt_31_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(srbd_wt_31_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")

###srbd run 31 alpha
srbd_alpha31bc1percents<-merge(subset(srbdoverall, srbdoverall$run == "31" & srbdoverall$variant == "Alpha"), srbd_alpha31bc1percents, by = "sample_seq", all= F)
nrow(srbd_alpha31bc1percents) #153 
table(srbd_alpha31bc1percents$V3 == srbd_alpha31bc1percents$seq_chars)
srbd_alpha31bc1percents<-unique(srbd_alpha31bc1percents)

srbd_alpha31bc2percents<-merge(subset(srbdoverall, srbdoverall$run == "31" & srbdoverall$variant == "Alpha"), srbd_alpha31bc2percents, by = "sample_seq", all= F)
nrow(srbd_alpha31bc2percents) #195 
table(srbd_alpha31bc2percents$V3 == srbd_alpha31bc2percents$seq_chars)
srbd_alpha31bc2percents<-unique(srbd_alpha31bc2percents)

srbd_alpha_31_info<-merge(srbd_alpha31bc1percents, srbd_alpha31bc2percents, by = "sample_seq", all = F)
head(srbd_alpha_31_info,3)
table(srbd_alpha_31_info$sampleID.x == srbd_alpha_31_info$sampleID.y) #64 
table(srbd_alpha_31_info$seq_chars.x == srbd_alpha_31_info$seq_chars.y)
srbd_alpha_31_info<-srbd_alpha_31_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(srbd_alpha_31_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")


##srbd run delta
srbd_deltabc1percents<-merge(subset(srbdoverall, srbdoverall$run == "139-141" & srbdoverall$variant == "Delta"), srbd_deltabc1percents, by = "sample_seq", all= F)
nrow(srbd_deltabc1percents) #125
table(srbd_deltabc1percents$V3 == srbd_deltabc1percents$seq_chars)
srbd_deltabc1percents<-unique(srbd_deltabc1percents)

srbd_deltabc2percents<-merge(subset(srbdoverall, srbdoverall$run == "139-141" & srbdoverall$variant == "Delta"), srbd_deltabc2percents, by = "sample_seq", all= F)
nrow(srbd_deltabc2percents) # 49
table(srbd_deltabc2percents$V3 == srbd_deltabc2percents$seq_chars)
srbd_deltabc2percents<-unique(srbd_deltabc2percents)

srbd_delta_info<-merge(srbd_deltabc1percents, srbd_deltabc2percents, by = "sample_seq", all = F)
head(srbd_delta_info,3)
table(srbd_delta_info$sampleID.x == srbd_delta_info$sampleID.y)
table(srbd_delta_info$seq_chars.x == srbd_delta_info$seq_chars.y) #15 now
srbd_delta_info<-srbd_delta_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(srbd_delta_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")


###srbd run omicron
srbd_omicronbc1percents<-merge(subset(srbdoverall, srbdoverall$run == "139-141" & srbdoverall$variant == "Omicron"), srbd_omicronbc1percents, by = "sample_seq", all= F)
nrow(srbd_omicronbc1percents) #573
table(srbd_omicronbc1percents$V3 == srbd_omicronbc1percents$seq_chars)
srbd_omicronbc1percents<-unique(srbd_omicronbc1percents)

srbd_omicronbc2percents<-merge(subset(srbdoverall, srbdoverall$run == "139-141" & srbdoverall$variant == "Omicron"), srbd_omicronbc2percents, by = "sample_seq", all= F)
nrow(srbd_omicronbc2percents) #441
table(srbd_omicronbc2percents$V3 == srbd_omicronbc2percents$seq_chars)
srbd_omicronbc2percents<-unique(srbd_omicronbc2percents)

srbd_omicron_info<-merge(srbd_omicronbc1percents, srbd_omicronbc2percents, by = "sample_seq", all = F)
head(srbd_omicron_info,3)
table(srbd_omicron_info$sampleID.x == srbd_omicron_info$sampleID.y) #62
table(srbd_omicron_info$seq_chars.x == srbd_omicron_info$seq_chars.y)
srbd_omicron_info<-srbd_omicron_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(srbd_omicron_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")

###all srbd
allsrbd<-rbind(srbd_wt_27_info, srbd_wt_31_info, srbd_alpha_31_info, srbd_delta_info, srbd_omicron_info)
nrow(allsrbd) #366 
length(unique(allsrbd$sampleID)) #231 


###combine all spbs too
#skip processing spbs 27 as we dropped the sequence from that output for being an artefact

###spbs 31
spbs_wt31bc1percents<-merge(subset(spbsoverall, spbsoverall$run == "31" & spbsoverall$variant == "WT"), spbs_wt31bc1percents, by = "sample_seq", all= F)
nrow(spbs_wt31bc1percents) #1
table(spbs_wt31bc1percents$V3 == spbs_wt31bc1percents$seq_chars)
spbs_wt31bc1percents<-unique(spbs_wt31bc1percents)

spbs_wt31bc2percents<-merge(subset(spbsoverall, spbsoverall$run == "31" & spbsoverall$variant == "WT"), spbs_wt31bc2percents, by = "sample_seq", all= F)
nrow(spbs_wt31bc2percents) #1
table(spbs_wt31bc2percents$V3 == spbs_wt31bc2percents$seq_chars)
spbs_wt31bc2percents<-unique(spbs_wt31bc2percents)

spbs_wt_31_info<-merge(spbs_wt31bc1percents, spbs_wt31bc2percents, by = "sample_seq", all = F)
head(spbs_wt_31_info,3)
table(spbs_wt_31_info$sampleID.x == spbs_wt_31_info$sampleID.y)
table(spbs_wt_31_info$seq_chars.x == spbs_wt_31_info$seq_chars.y)
spbs_wt_31_info<-spbs_wt_31_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(spbs_wt_31_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")

###spbs run 31 alpha
spbs_alpha31bc1percents<-merge(subset(spbsoverall, spbsoverall$run == "31" & spbsoverall$variant == "Alpha"), spbs_alpha31bc1percents, by = "sample_seq", all= F)
nrow(spbs_alpha31bc1percents) #4
table(spbs_alpha31bc1percents$V3 == spbs_alpha31bc1percents$seq_chars)
spbs_alpha31bc1percents<-unique(spbs_alpha31bc1percents)

spbs_alpha31bc2percents<-merge(subset(spbsoverall, spbsoverall$run == "31" & spbsoverall$variant == "Alpha"), spbs_alpha31bc2percents, by = "sample_seq", all= F)
nrow(spbs_alpha31bc2percents) #4
table(spbs_alpha31bc2percents$V3 == spbs_alpha31bc2percents$seq_chars)
spbs_alpha31bc2percents<-unique(spbs_alpha31bc2percents)

spbs_alpha_31_info<-merge(spbs_alpha31bc1percents, spbs_alpha31bc2percents, by = "sample_seq", all = F)
head(spbs_alpha_31_info,3)
table(spbs_alpha_31_info$sampleID.x == spbs_alpha_31_info$sampleID.y)
table(spbs_alpha_31_info$seq_chars.x == spbs_alpha_31_info$seq_chars.y)
spbs_alpha_31_info<-spbs_alpha_31_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(spbs_alpha_31_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")


##spbs run delta
spbs_deltabc1percents<-merge(subset(spbsoverall, spbsoverall$run == "139-141" & spbsoverall$variant == "Delta"), spbs_deltabc1percents, by = "sample_seq", all= F)
nrow(spbs_deltabc1percents) #7
table(spbs_deltabc1percents$V3 == spbs_deltabc1percents$seq_chars)
spbs_deltabc1percents<-unique(spbs_deltabc1percents)

spbs_deltabc2percents<-merge(subset(spbsoverall, spbsoverall$run == "139-141" & spbsoverall$variant == "Delta"), spbs_deltabc2percents, by = "sample_seq", all= F)
nrow(spbs_deltabc2percents) #7
table(spbs_deltabc2percents$V3 == spbs_deltabc2percents$seq_chars)
spbs_deltabc2percents<-unique(spbs_deltabc2percents)

spbs_delta_info<-merge(spbs_deltabc1percents, spbs_deltabc2percents, by = "sample_seq", all = F)
head(spbs_delta_info,3)
table(spbs_delta_info$sampleID.x == spbs_delta_info$sampleID.y)
table(spbs_delta_info$seq_chars.x == spbs_delta_info$seq_chars.y)
spbs_delta_info<-spbs_delta_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(spbs_delta_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")


###spbs run omicron
spbs_omicronbc1percents<-merge(subset(spbsoverall, spbsoverall$run == "139-141" & spbsoverall$variant == "Omicron"), spbs_omicronbc1percents, by = "sample_seq", all= F)
nrow(spbs_omicronbc1percents) #10
table(spbs_omicronbc1percents$V3 == spbs_omicronbc1percents$seq_chars)
spbs_omicronbc1percents<-unique(spbs_omicronbc1percents)

spbs_omicronbc2percents<-merge(subset(spbsoverall, spbsoverall$run == "139-141" & spbsoverall$variant == "Omicron"), spbs_omicronbc2percents, by = "sample_seq", all= F)
nrow(spbs_omicronbc2percents) #10
table(spbs_omicronbc2percents$V3 == spbs_omicronbc2percents$seq_chars)
spbs_omicronbc2percents<-unique(spbs_omicronbc2percents)

spbs_omicron_info<-merge(spbs_omicronbc1percents, spbs_omicronbc2percents, by = "sample_seq", all = F)
head(spbs_omicron_info,3)
table(spbs_omicron_info$sampleID.x == spbs_omicron_info$sampleID.y)
table(spbs_omicron_info$seq_chars.x == spbs_omicron_info$seq_chars.y)
spbs_omicron_info<-spbs_omicron_info[,c("sample_seq", "run.x", "variant.x", "sampleID.x", "seq_chars.x", "V4.x", "V4.y")]
colnames(spbs_omicron_info)<-c("sample_seq", "run", "variant", "sampleID", "seq_chars", "bc1_percent", "bc2_percent")


###all spbs
allspbs<-rbind(spbs_wt_31_info, spbs_alpha_31_info, spbs_delta_info, spbs_omicron_info)
nrow(allspbs) #22
length(unique(allspbs$sampleID)) #22


####now combine srbd and spbs
allsrbd$amplicon<-"S-Rbm"
allspbs$amplicon<-"S-Pbs"


###here we sum the percent per sample 
#we check if any samples are producing pQS at too high of a rate - any sample where over 5% of it is pQS is removed from the downstream analysis 
#as described in methods of manuscript
subset(allsrbd, allsrbd$bc1_percent > 5 | allsrbd$bc2_percent > 5)
#these are sequences where an individual sequence is over 5% of the sample. 
#                                                                                                   sample_seq     run variant  sampleID
#254 W62005079_ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTATTTTCCTTTACAATCATATGGTTTCCAACCCACTTATGGTGTT      31   Alpha W62005079
#310 X61404962_ATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT 139-141 Omicron X61404962
#314 X61404962_ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT 139-141 Omicron X61404962
#317 X61404962_ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT 139-141 Omicron X61404962
#355 X61701424_ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT 139-141 Omicron X61701424
#358 X61701424_ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT 139-141 Omicron X61701424
#                                                                                           seq_chars bc1_percent bc2_percent amplicon
#254 ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTATTTTCCTTTACAATCATATGGTTTCCAACCCACTTATGGTGTT   14.171254   10.170727    S-Rbm
#310 ATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT    4.635164    5.682795    S-Rbm
#314 ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT   16.622121   12.745319    S-Rbm
#317 ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT    4.943511    5.895033    S-Rbm
#355 ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT   25.254347   30.051646    S-Rbm
#358 ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT    4.765636    5.246787    S-Rbm

allsrbd[allsrbd$sampleID %in% c("W62005079", "X61404962", "X61701424" ),]

#check total percents per sample
allsrbd_summedpercent<-allsrbd %>% group_by(sampleID) %>% 
  summarise(sumBC1 = sum(bc1_percent), sumBC2 = sum(bc2_percent))
allsrbd_summedpercent[allsrbd_summedpercent$sumBC1 > 5 | allsrbd_summedpercent$sumBC2 > 5,]
# sampleID  sumBC1 sumBC2
#  <chr>      <dbl>  <dbl>
#1 W62005079   14.2   10.2
#2 X61404962   35.2   37.4
#3 X61701424   43.6   49.9
#still only these three samples.


###we drop the 3 samples with too high percents##
allsrbd_exclude3<-subset(allsrbd, !(allsrbd$sampleID %in% c("W62005079", "X61404962", "X61701424")))
subset(allsrbd_exclude3, allsrbd_exclude3$bc1_percent > 5 | allsrbd_exclude3$bc2_percent > 5)

allsrbdspbs<-rbind(allsrbd_exclude3, allspbs)
dim(allsrbdspbs) #363 8

#save this table
write.xlsx(allsrbdspbs, file= "aggregated_QS_percents_SrbmAndSpbs_dropped3Samples_ForSuppFigs.xlsx")
allsrbdspbs<-read.xlsx("aggregated_QS_percents_SrbmAndSpbs_dropped3Samples_ForSuppFigs.xlsx")

mean(allsrbdspbs$bc1_percent) #0.3764408
mean(allsrbdspbs$bc2_percent) #0.3108724

#sd
sd(allsrbdspbs$bc1_percent) #0.6763237
sd(allsrbdspbs$bc2_percent) #0.2959845
#reported these numbers to manuscript


##allsrbdspbs
head(allsrbdspbs, 3)
allsrbdspbs_forplot<-allsrbdspbs[,c("variant",  "sampleID", "bc1_percent", "bc2_percent", "amplicon")]
allsrbdspbs_forplotwide<-melt(allsrbdspbs_forplot, id.vars=c("variant", "amplicon", "sampleID"))
allsrbdspbs_forplotwide$BC<-"Barcode Copy 1"
allsrbdspbs_forplotwide[allsrbdspbs_forplotwide$variable == "bc2_percent",]$BC<-"Barcode Copy 2"
table(allsrbdspbs_forplotwide$BC)
allsrbdspbs_forplotwide$variable<-NULL
head(allsrbdspbs_forplotwide)
allsrbdspbs_forplotwide$variant_F<-factor(allsrbdspbs_forplotwide$variant, levels=c("WT", "Alpha", "Delta", "Omicron"))
allsrbdspbs_forplotwide$amplicon_F<-factor(allsrbdspbs_forplotwide$amplicon, levels=c("S-Rbm","S-Pbs"))


###next make table of num samples, num seqs, num samples that produced QS, rate of qs production
nrow(allsrbd_exclude3) #341 
nrow(subset(allsrbd_exclude3, allsrbd_exclude3$bc1_percent < 5.00 & allsrbd_exclude3$bc2_percent < 5.00))#341, all are now under 5%
length(unique(allsrbd_exclude3$sampleID)) #228 samples,
length(unique(allsrbd_exclude3$sample_seq)) #341

nrow(allspbs) #22 
nrow(subset((allspbs), (allspbs)$bc1_percent < 5.00 & (allspbs)$bc2_percent < 5.00)) #22 out of 22 so 100%
##
length(unique(allspbs$sampleID)) #21 samples,
length(unique(allspbs$sample_seq)) #22


###save tables for downstream analysis
write.xlsx(allsrbd_exclude3[,c("run",	"variant",	"sampleID",	"seq_chars")], file = "overall_srbd_table.xlsx")
write.xlsx(allspbs[,c("run",	"variant",	"sampleID",	"seq_chars")], file = "overall_spbs_table.xlsx")

######also make boxplot of total qs seq per sample ####
#so like if sample 1 had qs.1 with 3 % and qs.2 with 3 % it has 6% total 
allsrbdspbs<-read.xlsx("aggregated_QS_percents_SrbmAndSpbs_dropped3Samples_ForSuppFigs.xlsx")
head(allsrbdspbs)
allsrbdspbs_sumpersample<-as.data.frame(allsrbdspbs %>% group_by(variant, sampleID, amplicon) %>% summarise(totalBC1percent = sum(bc1_percent), totalBC2percent = sum(bc2_percent)))


allsrbdspbs_sumpersamplewide<-melt(allsrbdspbs_sumpersample, id.vars=c("variant", "amplicon", "sampleID"))
allsrbdspbs_sumpersamplewide$BC<-"Barcode Copy 1"
allsrbdspbs_sumpersamplewide[allsrbdspbs_sumpersamplewide$variable == "totalBC2percent",]$BC<-"Barcode Copy 2"
head(allsrbdspbs_sumpersamplewide)

allsrbdspbs_sumpersamplewide$variant_F<-factor(allsrbdspbs_sumpersamplewide$variant, levels=c("WT", "Alpha", "Delta", "Omicron"))
allsrbdspbs_sumpersamplewide$amplicon_F<-factor(allsrbdspbs_sumpersamplewide$amplicon, levels=c("S-Rbm","S-Pbs"))

##plot this
bc_boxplot<-ggplot(allsrbdspbs_sumpersamplewide, aes(x=variant_F, y=value, colour = BC)) + facet_wrap(.~amplicon_F) +
    geom_jitter(size = 0.5, alpha= 0.5, position = position_jitterdodge()) + ggtitle("Quasispecies: Total Percent per Sample\nthat is QS sequences") + 
  xlab("Variant") + ylab("Percent") + theme_bw() + ylim(0,5) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(size = 7,  hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 8))
bc_boxplot

pdf("SuppFig13_Boxplot_QSPercentSumPerSample.pdf", width = 7, height = 5)
bc_boxplot
dev.off()


###then remake the table of sample - seq left after filtering, then the bubble plots and the final sequence alignments 
srbduniqsample<-table(allsrbd_exclude3$sampleID)
length(srbduniqsample[srbduniqsample > 1]) #109

table(allsrbd_exclude3$variant) #count of seqs
# Alpha   Delta Omicron      WT 
#     63      15      38     225 >> after removing 3 high % samples
     
table(allspbs$variant)
#Alpha   Delta Omicron      WT 
#      4       7      10       1 

length(unique(allsrbd_exclude3$sampleID)) #228
length(unique(allsrbd_exclude3[allsrbd_exclude3$variant == "WT",]$sampleID)) # 145
length(unique(allsrbd_exclude3[allsrbd_exclude3$variant == "Alpha",]$sampleID)) # 38
length(unique(allsrbd_exclude3[allsrbd_exclude3$variant == "Delta",]$sampleID)) # 12
length(unique(allsrbd_exclude3[allsrbd_exclude3$variant == "Omicron",]$sampleID)) # 33

###srbm samples with multiple seqs
srbm_multiseq<-table(allsrbd_exclude3$sampleID)
srbm_multiseq<-srbm_multiseq[order(srbm_multiseq)] #
length(srbm_multiseq[srbm_multiseq > 1]) #109 samples with multiple seq 


##unique srbm sequences per var after excluding 3 samples 
length(table(allsrbd_exclude3[allsrbd_exclude3$variant == "WT",]$seq_chars)) #13
length(table(allsrbd_exclude3[allsrbd_exclude3$variant == "Alpha",]$seq_chars)) #5
length(table(allsrbd_exclude3[allsrbd_exclude3$variant == "Delta",]$seq_chars)) #6
length(table(allsrbd_exclude3[allsrbd_exclude3$variant == "Omicron",]$seq_chars)) #5

###spbs samples with multiple seqs
spbs_multiseq<-table(allspbs$sampleID)
spbs_multiseq<-spbs_multiseq[order(spbs_multiseq)] #
length(spbs_multiseq[spbs_multiseq > 1]) #1 samples with multiple seq 


length(unique(allspbs$sampleID)) #21
length(unique(allspbs[allspbs$variant == "WT",]$sampleID)) #1
length(unique(allspbs[allspbs$variant == "Alpha",]$sampleID)) #4
length(unique(allspbs[allspbs$variant == "Delta",]$sampleID)) #6
length(unique(allspbs[allspbs$variant == "Omicron",]$sampleID)) #10

table(unique(allsrbd_exclude3$sampleID) %in% unique(allspbs$sampleID) ) #unique: 219 not in both, 9 in both

###remake bubble plots#### 
library(Biostrings)
#output list of kept seqs to edit the other inputs
allsrbd_exclude3_variantSeqs<-allsrbd_exclude3[,c("variant", "seq_chars")]
allsrbd_exclude3_variantSeqs<-unique(allsrbd_exclude3_variantSeqs)
common_dnastring1<-DNAStringSet(allsrbd_exclude3_variantSeqs$seq_chars)
common_aas1<-as.data.frame(translate(common_dnastring1))
allsrbd_exclude3_variantSeqs<-cbind(allsrbd_exclude3_variantSeqs, common_aas1)
colnames(allsrbd_exclude3_variantSeqs)<-c("variant", "seq_chars",  "AAstring")
write.xlsx(allsrbd_exclude3_variantSeqs, file = "srbd_unique_kept_sequences_for_figs.xlsx")
allsrbd_exclude3_variantSeqs<-read.xlsx("srbd_unique_kept_sequences_for_figs.xlsx")

allspbs_variantSeqs<-allspbs[,c("variant", "seq_chars")]
allspbs_variantSeqs<-unique(allspbs_variantSeqs)
common_dnastring2<-DNAStringSet(allspbs_variantSeqs$seq_chars)
common_aas2<-as.data.frame(translate(common_dnastring2))
#can be a warning here for dropping the last 2 bases but that's fine
allspbs_variantSeqs<-cbind(allspbs_variantSeqs, common_aas2)
colnames(allspbs_variantSeqs)<-c("variant", "seq_chars",  "AAstring")
head(allspbs_variantSeqs)
write.xlsx(allspbs_variantSeqs, file = "spbs_unique_kept_sequences_for_figs.xlsx")
allspbs_variantSeqs<-read.xlsx("spbs_unique_kept_sequences_for_figs.xlsx")


#### make bar graphs####
allsrbd_exclude3
length(unique(allsrbd_exclude3$sampleID)) #228
length(unique( allspbs$sampleID)) #21


##total number samples = 
length(unique(c(allsrbd_exclude3$sampleID, allspbs$sampleID)))
#240 unique samples 
length(unique(c(allsrbd_exclude3$seq_chars, allspbs$seq_chars)))
##34 unique sequences
nrow(allsrbd_exclude3) #341
nrow(allspbs) #22
##so total qs seq is 341+ 22 = 363 and total unique samples involved is 240 

allsrbd_exclude3_grouped<-as.data.frame(allsrbd_exclude3 %>% group_by(amplicon, variant) %>% summarise(countofseq = n()))
#sum is 225+38+15+63 so 341
allspbs_grouped<-as.data.frame(allspbs %>% group_by(amplicon, variant) %>% summarise(countofseq = n()))
#sum is 22
#so 341+22 = 363 

allsrbd_exclude3_grouped$percent<-allsrbd_exclude3_grouped$countofseq / 363*100
allsrbd_exclude3_grouped$variant_F<-factor(allsrbd_exclude3_grouped$variant, levels=c("WT", "Alpha", "Delta", "Omicron" ))

allspbs_grouped$percent<-allspbs_grouped$countofseq / 363*100
allspbs_grouped$variant_F<-factor(allspbs_grouped$variant, levels=c("WT", "Alpha", "Delta", "Omicron" ))

#make bar graph
srbd_exclude3_bars<-ggplot(allsrbd_exclude3_grouped, aes(variant_F, percent, fill = variant)) + geom_col(colour = "black") +
  ggtitle("S-Rbm pQS (341/363), percent of all pQS by VoC\n[total pQS = 363 sequences found in 240 samples]") +  ylab("Percent") + scale_x_discrete(drop=FALSE) +
  xlab("Week") + theme_bw() + scale_y_continuous(limits=c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),   axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 14))
spbs_bars<-ggplot(allspbs_grouped, aes(variant_F, percent, fill = variant)) + geom_col(colour = "black") +
  ggtitle("S-Pbs pQS (22/363), percent of all pQS by VoC\n[total pQS = 363 sequences found in 240 samples]") +  ylab("Percent") + scale_x_discrete(drop=FALSE) +
  xlab("Week") + theme_bw() + scale_y_continuous(limits=c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),   axis.text.x = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 14))

pdf("qsBars_Fig6_andSup12.pdf", height = 5, width = 6)
srbd_exclude3_bars
spbs_bars
dev.off()


###bubble plot code####
#hardcode references
srbd_omicron<-"ATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT"#GGTC"
srbd_omicron_AA<-"IYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGV"#G"
srbd_delta<-"ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"#GGTT"
srbd_delta_AA<-"IYQAGSKPCNGVEGFNCYFPLQSYGFQPTNGV"#G"
srbd_wt<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"#GGTT"
srbd_wt_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGV"#G"
srbd_alpha<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTTATGGTGTT"
srbd_alpha_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGV"


#function to splice up sequences
splitpos_control<-function(x){
  outputdf <- data.frame(position=character(), base=character(), stringsAsFactors=FALSE) 
  colnames(outputdf)<-c("position", "base")
  for(j in 1:nchar(x)){
    base<-substring(x, j, j)
    outputdf[nrow(outputdf)+1,]<-c( j, base)
  }
  return(outputdf)
}


srbdwtaa_split<-splitpos_control(srbd_wt_AA)
srbdalphaaa_split<-splitpos_control(srbd_alpha_AA)
srbddeltaaa_split<-splitpos_control(srbd_delta_AA)
srbdomicronaa_split<-splitpos_control(srbd_omicron_AA)

colnames(srbdwtaa_split)<-c("position", "refseq_AA")
colnames(srbdalphaaa_split)<-c("position", "refseq_AA")
colnames(srbddeltaaa_split)<-c("position", "refseq_AA")
colnames(srbdomicronaa_split)<-c("position", "refseq_AA")


#also our layout file to connect info and ordering to the AA seq 
#we made this file manually using our pQS outputs to be sure our sequences were labeled clearly
layouts<-read.xlsx("srbd_qs_aa_id_layout.xlsx")
order_of_seqs<-layouts$seqname

split_pos_layout <- function(x){
  outputdf <- data.frame(run_variant=character(),  annotation=character(), seqname=character(),
                 position=character(),  base=character(), stringsAsFactors=FALSE) 
  colnames(outputdf)<-c("run_variant", "annotation", "seqname", "position", "base")
  for(i in 1:nrow(x)){ #for each row of table
    print(i)
    annotation<-x[i,"annotation" ]
    run_variant<-x[i,"run_variant" ]
    seqname<-x[i,"seqname"]
    aasequence<-x[i,"aasequence" ]
    #
    for(j in 1:nchar(aasequence)){
      base<-substring(aasequence, j, j)
      outputdf[nrow(outputdf)+1,]<-c(run_variant, annotation, seqname, j, base)
      
    }
  }
  return(outputdf)
}

processedlayout<-split_pos_layout(layouts)

#split into the 4 types and merge to correct refseq bases
#then mark change or no change
processedlayout_wt<-subset(processedlayout, processedlayout$run_variant == "wt")
processedlayout_alpha<-subset(processedlayout, processedlayout$run_variant == "alpha")
processedlayout_delta<-subset(processedlayout, processedlayout$run_variant == "delta")
processedlayout_omicron<-subset(processedlayout, processedlayout$run_variant == "omicron")

processedlayout_wt<-merge(processedlayout_wt, srbdwtaa_split, by = "position" )
processedlayout_alpha<-merge(processedlayout_alpha, srbdalphaaa_split, by = "position" )
processedlayout_delta<-merge(processedlayout_delta, srbddeltaaa_split, by = "position" )
processedlayout_omicron<-merge(processedlayout_omicron, srbdomicronaa_split, by = "position" )

all_processed<-rbind(processedlayout_wt, processedlayout_alpha, processedlayout_delta, processedlayout_omicron)
all_processed$change<-"none"

all_processed[all_processed$refseq_AA != all_processed$base,]$change<-all_processed[all_processed$refseq_AA != all_processed$base,]$base
all_processed$magnitude<-stringr::str_extract(all_processed$annotation, "\\([0-9]{1,3}\\)")
all_processed$magnitude<-stringr::str_extract(all_processed$magnitude, "[0-9]{1,3}")
all_processed_changeonly<-subset(all_processed, all_processed$change != "none")

#manage Synonymous sequences - since they have no change to remain in the list 
all_processed_synonymous<-subset(all_processed, all_processed$seqname %in% c(
  "WT_1_Synonymous",  "Alpha_1_Synonymous", "Delta_1_Synonymous.1",  "Delta_1_Synonymous.2", "Delta_1_Synonymous.3", "Omicron_15_Synonymous" ))
#set magnitudes for 0
all_processed_synonymous$magnitude<-"NA"

all_processed_changeonly<-rbind(all_processed_changeonly, all_processed_synonymous)
all_processed_changeonly$magnitude<-as.numeric(all_processed_changeonly$magnitude)
all_processed_changeonly$position<-as.numeric(all_processed_changeonly$position)
all_processed_changeonly$annotation_F<-factor(all_processed_changeonly$seqname, levels=rev(order_of_seqs))
all_processed_changeonly$run_variant<-factor(all_processed_changeonly$run_variant, levels = c("wt", "alpha", "delta", "omicron"))#levels = c("omicron", "delta", "alpha", "wt"))

#fix x axis numbers $position+471
all_processed_changeonlyAA_v2<-all_processed_changeonly
all_processed_changeonlyAA_v2$position<-all_processed_changeonlyAA_v2$position+471

all_processed_changeonlyAA_v2$magnitude
all_processed_changeonlyAA_v3<-all_processed_changeonlyAA_v2
all_processed_changeonlyAA_v3[!is.na(all_processed_changeonlyAA_v3$magnitude) & all_processed_changeonlyAA_v3$magnitude > 50,]$magnitude<-50


pdf("srbd_qs_dotplot_AA.pdf", width = 10, height = 9)
aapanel<-ggplot(all_processed_changeonlyAA_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y") +
  xlab("SRBD SEQUENCE POSITION") + ylab("Sequence ID") + theme_bw() + ggtitle("SRBD Sequence Positions - AA QS Changes") +
  scale_x_continuous(limits = c(470,505), breaks = seq(470,505,5)) + scale_size_area(breaks = c(1, 5, 10, 25, 50), limits = c(0, 50)) +
  labs(subtitle = "Individual Refseqs Used", fill = "AA Change") +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 12))
aapanel
dev.off()

aapanel_no_yaxislab<-ggplot(all_processed_changeonlyAA_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() +  facet_grid(run_variant~., scales = "free_y", space = "free_y") +
  xlab("SRBD AA SEQUENCE POSITION") +  theme_bw() +  scale_x_continuous(limits = c(470,505), breaks = seq(470,505,5)) +
  scale_size_area(breaks = c(1, 5, 10, 25, 50), limits = c(0, 50)) + labs(fill = "AA Change") +
  theme(axis.title.x = element_text(size=12), axis.title.y = element_blank(), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, hjust = 0.5), axis.ticks.y = element_blank(),
      axis.text.y = element_blank())
aapanel_no_yaxislab


####also repeat for NTDs so we can put them side by side ########
srbdwtsplit<-splitpos_control(srbd_wt)
srbdalphasplit<-splitpos_control(srbd_alpha)
srbddeltasplit<-splitpos_control(srbd_delta)
srbdomicronsplit<-splitpos_control(srbd_omicron)

colnames(srbdwtsplit)<-c("position", "refseq")
colnames(srbdalphasplit)<-c("position", "refseq")
colnames(srbddeltasplit)<-c("position", "refseq")
colnames(srbdomicronsplit)<-c("position", "refseq")

#also our manually made layout file to connect info and ordering to the AA seq 
layouts<-read.xlsx("srbd_qs_aa_id_layout.xlsx")
order_of_seqs<-layouts$seqname

#crunch sequences
split_pos_layout <- function(x){
  outputdf <- data.frame(run_variant=character(),  annotation=character(), seqname=character(),
                 position=character(),  base=character(), stringsAsFactors=FALSE) 
  colnames(outputdf)<-c("run_variant", "annotation", "seqname", "position", "base")
  for(i in 1:nrow(x)){ #for each row of table
    print(i)
    annotation<-x[i,"annotation" ]
    run_variant<-x[i,"run_variant" ]
    seqname<-x[i,"seqname"]
    ntdsequence<-x[i,"ntdsequence" ]
    #
    for(j in 1:nchar(ntdsequence)){
      base<-substring(ntdsequence, j, j)
      outputdf[nrow(outputdf)+1,]<-c(run_variant, annotation, seqname, j, base)
      
    }
  }
  return(outputdf)
}

processedlayout<-split_pos_layout(layouts)

#split into the 4 types and merge to correct refseq bases
#then mark change or no change
processedlayout_wt<-subset(processedlayout, processedlayout$run_variant == "wt")
processedlayout_alpha<-subset(processedlayout, processedlayout$run_variant == "alpha")
processedlayout_delta<-subset(processedlayout, processedlayout$run_variant == "delta")
processedlayout_omicron<-subset(processedlayout, processedlayout$run_variant == "omicron")

processedlayout_wt<-merge(processedlayout_wt, srbdwtsplit, by = "position" )
processedlayout_alpha<-merge(processedlayout_alpha, srbdalphasplit, by = "position" )
processedlayout_delta<-merge(processedlayout_delta, srbddeltasplit, by = "position" )
processedlayout_omicron<-merge(processedlayout_omicron, srbdomicronsplit, by = "position" )

all_processed<-rbind(processedlayout_wt, processedlayout_alpha, processedlayout_delta, processedlayout_omicron)
all_processed$change<-"none"
all_processed[all_processed$refseq != all_processed$base,]$change<-all_processed[all_processed$refseq != all_processed$base,]$base
all_processed$magnitude<-stringr::str_extract(all_processed$annotation, "\\([0-9]{1,3}\\)")
all_processed$magnitude<-stringr::str_extract(all_processed$magnitude, "[0-9]{1,3}")


#adjust synonymous
all_processed_changeonly<-subset(all_processed, all_processed$change != "none")
all_processed_synonymous<-subset(all_processed, all_processed$seqname %in% c(
   "WT_1_Synonymous", "Alpha_1_Synonymous", "Delta_1_Synonymous.1",
  "Delta_1_Synonymous.2", "Delta_1_Synonymous.3" ))
all_processed_synonymous$magnitude<-0
all_processed_synonymous$magnitude<-"NA"

all_processed_changeonly<-rbind(all_processed_changeonly, all_processed_synonymous)
all_processed_changeonly$magnitude<-as.numeric(all_processed_changeonly$magnitude)
all_processed_changeonly$position<-as.numeric(all_processed_changeonly$position)
all_processed_changeonly$annotation_F<-factor(all_processed_changeonly$seqname, levels=rev(order_of_seqs))
all_processed_changeonly$run_variant<-factor(all_processed_changeonly$run_variant, levels = c("wt", "alpha", "delta", "omicron"))#levels = c("omicron", "delta", "alpha", "wt"))

#adjust x axis numbers $position+471
all_processed_changeonly_v2<-all_processed_changeonly
all_processed_changeonly_v2$position<-all_processed_changeonly_v2$position+1413
all_processed_changeonly_v3<-all_processed_changeonly_v2
all_processed_changeonly_v3[!is.na(all_processed_changeonly_v3$magnitude) & all_processed_changeonly_v3$magnitude > 50,]$magnitude<-50

pdf("srbd_qs_dotplot_NTD.pdf", width = 10, height = 9)
ntdpanel<-ggplot(all_processed_changeonly_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y", switch = "y") +
  xlab("SRBD SEQUENCE POSITION") + ylab("Sequence ID") + theme_bw() + ggtitle("SRBD NTD QS Changes") + 
  scale_x_continuous(limits = c(1413,1510), breaks = seq(1410,1510,10)) +  scale_size_area(breaks = c(1, 5, 10, 25, 50), limits = c(0, 50)) +
  labs(subtitle = "Individual Refseqs Used", fill = "NTD Change") +
  theme(axis.title = element_text(size=12), legend.position = c("left"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 10))
ntdpanel
dev.off()

ntdpanel_forpatch<-ggplot(all_processed_changeonly_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) +
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y", switch = "y") +  xlab("SRBD NTD SEQUENCE POSITION") +
  ylab("Sequence ID") + theme_bw() +  scale_x_continuous(limits = c(1413,1510), breaks = seq(1410,1510,10)) +
  scale_size_area(breaks = c(1, 5, 10, 25, 50), limits = c(0, 50)) +  labs(fill = "NTD Change") +
  theme(axis.title = element_text(size=12), legend.position = c("left"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8, hjust = 0.5), axis.text.y = element_text(hjust = 1, size = 10))
ntdpanel_forpatch


###combo panel####
library(patchwork)
patch <- ntdpanel_forpatch | aapanel_no_yaxislab 
patch<- patch + plot_annotation(title = 'SRBD QS Changes (NTDs on left and AAs on right)', theme = theme(plot.title = element_text(hjust = 0.5)))

pdf("combo_bubbles_SRBD_NTDandAA.pdf", width = 14, height = 12)
patch
dev.off()



####SPBS bubbles
spbs_wt<-"AATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACT"#GGTC"
spbs_wt_AA<-"NSPRRARSVASQSIIAYTMS"#G"
spbs_alpha<-"AATTCTCATCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACT"#GGTT"
spbs_alpha_AA<-"NSHRRARSVASQSIIAYTMS"#G"
spbs_delta<-"AATTCTCGTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACT"#GGTT"
spbs_delta_AA<-"NSRRRARSVASQSIIAYTMS"#G"
spbs_omicron<-"AAGTCTCATCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACT"
spbs_omicron_AA<-"KSHRRARSVASQSIIAYTMS"


#wt seq dataframe
splitpos_control<-function(x){
  outputdf <- data.frame(position=character(), base=character(), stringsAsFactors=FALSE) 
  colnames(outputdf)<-c("position", "base")
  for(j in 1:nchar(x)){
    base<-substring(x, j, j)
    outputdf[nrow(outputdf)+1,]<-c( j, base)
  }
  return(outputdf)
}

spbswtaa_split<-splitpos_control(spbs_wt_AA)
spbsalphaaa_split<-splitpos_control(spbs_alpha_AA)
spbsdeltaaa_split<-splitpos_control(spbs_delta_AA)
spbsomicronaa_split<-splitpos_control(spbs_omicron_AA)

colnames(spbswtaa_split)<-c("position", "refseq_AA")
colnames(spbsalphaaa_split)<-c("position", "refseq_AA")
colnames(spbsdeltaaa_split)<-c("position", "refseq_AA")
colnames(spbsomicronaa_split)<-c("position", "refseq_AA")

##import our data from spbs analysis
alldata<-read.xlsx("overall_spbs_table.xlsx")
head(alldata)

#also our layout file to connect info and ordering to the AA seq - made this file manually to be consistent with naming
layouts<-read.xlsx("spbs_qs_aa_id_layout.xlsx")
unique(layouts$aasequence)
#there are dupes in the sequences due to multiple synonymous seqs and then due to some reversions matching other var seqs
order_of_seqs<-layouts$seqname

#crunch positions
split_pos_layout <- function(x){
  outputdf <- data.frame(run_variant=character(),  annotation=character(), seqname=character(),
                 position=character(),  base=character(), stringsAsFactors=FALSE) 
  colnames(outputdf)<-c("run_variant", "annotation", "seqname", "position", "base")
  for(i in 1:nrow(x)){ #for each row of table
    print(i)
    annotation<-x[i,"annotation" ]
    run_variant<-x[i,"run_variant" ]
    seqname<-x[i,"seqname"]
    aasequence<-x[i,"aasequence" ]
    #
    for(j in 1:nchar(aasequence)){
      base<-substring(aasequence, j, j)
      outputdf[nrow(outputdf)+1,]<-c(run_variant, annotation, seqname, j, base)
      
    }
  }
  return(outputdf)
}

processedlayout<-split_pos_layout(layouts)

#split into the 4 types and merge to correct refseq bases
#then mark change or no change
processedlayout_wt<-subset(processedlayout, processedlayout$run_variant == "wt")
processedlayout_alpha<-subset(processedlayout, processedlayout$run_variant == "alpha")
processedlayout_delta<-subset(processedlayout, processedlayout$run_variant == "delta")
processedlayout_omicron<-subset(processedlayout, processedlayout$run_variant == "omicron")

processedlayout_wt<-merge(processedlayout_wt, spbswtaa_split, by = "position" )
processedlayout_alpha<-merge(processedlayout_alpha, spbsalphaaa_split, by = "position" )
processedlayout_delta<-merge(processedlayout_delta, spbsdeltaaa_split, by = "position" )
processedlayout_omicron<-merge(processedlayout_omicron, spbsomicronaa_split, by = "position" )

all_processed<-rbind(processedlayout_wt, processedlayout_alpha, processedlayout_delta, processedlayout_omicron)
all_processed$change<-"none"
all_processed[all_processed$refseq_AA != all_processed$base,]$change<-all_processed[all_processed$refseq_AA != all_processed$base,]$base

all_processed$magnitude<-stringr::str_extract(all_processed$annotation, "\\([0-9]{1,3}\\)")
all_processed$magnitude<-stringr::str_extract(all_processed$magnitude, "[0-9]{1,3}")

all_processed_changeonly<-subset(all_processed, all_processed$change != "none")
#set up synonymous ones
all_processed_synonymous<-subset(all_processed, all_processed$seqname %in% c("Delta_2_Synonymous"))
all_processed_synonymous$magnitude<-"NA"

all_processed_changeonly<-rbind(all_processed_changeonly, all_processed_synonymous)
all_processed_changeonly$magnitude<-as.numeric(all_processed_changeonly$magnitude)
all_processed_changeonly$position<-as.numeric(all_processed_changeonly$position)
all_processed_changeonly$annotation_F<-factor(all_processed_changeonly$seqname, levels=rev(order_of_seqs))
all_processed_changeonly$run_variant<-factor(all_processed_changeonly$run_variant, levels = c("wt", "alpha", "delta", "omicron"))#levels = c("omicron", "delta", "alpha", "wt"))

#adjust x axis numbers $position+471
all_processed_changeonlyAA_v2<-all_processed_changeonly
all_processed_changeonlyAA_v2$position<-all_processed_changeonlyAA_v2$position+678 #first N in spbs should be 679 and should go to 98
all_processed_changeonlyAA_v3<-all_processed_changeonlyAA_v2


pdf("qs_dotplot_AA_spbs.pdf", width = 10, height = 8)
aapanel<-ggplot(all_processed_changeonlyAA_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y") +  xlab("SPBS SEQUENCE POSITION") + ylab("Sequence ID") + theme_bw() + 
  ggtitle("SPBS Sequence Positions - AA QS Changes") +  scale_x_continuous(limits = c(675,699), breaks = seq(675,699,5)) +
  scale_size_area(breaks = c(1, 5, 10, 15, 20), limits = c(0, 20)) +  labs(subtitle = "Individual Refseqs Used", fill = "AA Change") +
  theme(axis.title = element_text(size=12), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 12))
aapanel
dev.off()

aapanel_no_yaxislab<-ggplot(all_processed_changeonlyAA_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y") +  theme_bw() + ggtitle("SPBS Sequence Positions - AA QS Changes") + 
  scale_x_continuous(limits = c(675,699), breaks = seq(675,699,5)) +  scale_size_area(breaks = c(1, 5, 10, 15, 20), limits = c(0, 20)) +  labs(fill = "AA Change") +
  theme(axis.title.x = element_text(size=12), axis.title.y = element_blank(), legend.position = c("right"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, hjust = 0.5), axis.ticks.y = element_blank(),
      axis.text.y = element_blank())
aapanel_no_yaxislab


####also repeat for NTDs so we can put them side by side ########
spbswtsplit<-splitpos_control(spbs_wt)
spbsalphasplit<-splitpos_control(spbs_alpha)
spbsdeltasplit<-splitpos_control(spbs_delta)
spbsomicronsplit<-splitpos_control(spbs_omicron)

colnames(spbswtsplit)<-c("position", "refseq")
colnames(spbsalphasplit)<-c("position", "refseq")
colnames(spbsdeltasplit)<-c("position", "refseq")
colnames(spbsomicronsplit)<-c("position", "refseq")

##import our data from spbs analysis
alldata<-read.xlsx("overall_spbs_table.xlsx")
head(alldata)

#also our layout file to connect info and ordering to the AA seq - we made this manually just to have labels organized
layouts<-read.xlsx("spbs_qs_aa_id_layout.xlsx")
layouts
unique(layouts$ntdsequence)
#there are dupes in the sequences due to multiple synonymous seqs and then due to some reversions matching other var seqs
order_of_seqs<-layouts$seqname

#crunch positions
split_pos_layout <- function(x){
  outputdf <- data.frame(run_variant=character(),  annotation=character(), seqname=character(),
                 position=character(),  base=character(), stringsAsFactors=FALSE) 
  colnames(outputdf)<-c("run_variant", "annotation", "seqname", "position", "base")
  for(i in 1:nrow(x)){ #for each row of table
    print(i)
    annotation<-x[i,"annotation" ]
    run_variant<-x[i,"run_variant" ]
    seqname<-x[i,"seqname"]
    ntdsequence<-x[i,"ntdsequence" ]
    #
    for(j in 1:nchar(ntdsequence)){
      base<-substring(ntdsequence, j, j)
      outputdf[nrow(outputdf)+1,]<-c(run_variant, annotation, seqname, j, base)
      
    }
  }
  return(outputdf)
}

processedlayout<-split_pos_layout(layouts)

#split into the 4 types and merge to correct refseq bases
#then mark change or no change
processedlayout_wt<-subset(processedlayout, processedlayout$run_variant == "wt")
processedlayout_alpha<-subset(processedlayout, processedlayout$run_variant == "alpha")
processedlayout_delta<-subset(processedlayout, processedlayout$run_variant == "delta")
processedlayout_omicron<-subset(processedlayout, processedlayout$run_variant == "omicron")

processedlayout_wt<-merge(processedlayout_wt, spbswtsplit, by = "position" )
processedlayout_alpha<-merge(processedlayout_alpha, spbsalphasplit, by = "position" )
processedlayout_delta<-merge(processedlayout_delta, spbsdeltasplit, by = "position" )
processedlayout_omicron<-merge(processedlayout_omicron, spbsomicronsplit, by = "position" )

all_processed<-rbind(processedlayout_wt, processedlayout_alpha, processedlayout_delta, processedlayout_omicron)
all_processed$change<-"none"

all_processed[all_processed$refseq != all_processed$base,]$change<-all_processed[all_processed$refseq != all_processed$base,]$base
all_processed$magnitude<-stringr::str_extract(all_processed$annotation, "\\([0-9]{1,3}\\)")
all_processed$magnitude<-stringr::str_extract(all_processed$magnitude, "[0-9]{1,3}")


#need to set up the synonymous rows 
all_processed_changeonly<-subset(all_processed, all_processed$change != "none")
all_processed_synonymous<-subset(all_processed, all_processed$seqname %in% c( "Delta_2_Synonymous", "Omicron_20_Synonymous" ))
all_processed_synonymous$magnitude<-0
all_processed_synonymous$magnitude<-"NA"

all_processed_changeonly<-rbind(all_processed_changeonly, all_processed_synonymous)
all_processed_changeonly$magnitude<-as.numeric(all_processed_changeonly$magnitude)
all_processed_changeonly$position<-as.numeric(all_processed_changeonly$position)
all_processed_changeonly$annotation_F<-factor(all_processed_changeonly$seqname, levels=rev(order_of_seqs))
all_processed_changeonly$run_variant<-factor(all_processed_changeonly$run_variant, levels = c("wt", "alpha", "delta", "omicron"))

#replace axis numbers
all_processed_changeonly_v2<-all_processed_changeonly
all_processed_changeonly_v2$position<-all_processed_changeonly_v2$position+2034
all_processed_changeonly_v3<-all_processed_changeonly_v2

pdf("qs_dotplot_NTD_spbs.pdf", width = 10, height = 8)
ntdpanel<-ggplot(all_processed_changeonly_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y", switch = "y") + xlab("SPBS SEQUENCE POSITION") + ylab("Sequence ID") + theme_bw() + 
  ggtitle("SPBS NTD QS Changes") + scale_x_continuous(limits = c(2034,2100), breaks = seq(2035,2095,10)) +
  scale_size_area(breaks = c(1, 5, 10, 15, 20), limits = c(0, 20)) +  labs(subtitle = "Individual Refseqs Used", fill = "NTD Change") +
  theme(axis.title = element_text(size=12), legend.position = c("left"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 10))
ntdpanel
dev.off()

ntdpanel_forpatch<-ggplot(all_processed_changeonly_v3, aes(x=position, y=annotation_F, colour = base, size = magnitude)) + 
  geom_point() + facet_grid(run_variant~., scales = "free_y", space = "free_y", switch = "y") +
  xlab("SPBS NTD SEQUENCE POSITION") + ylab("Sequence ID") + theme_bw() + scale_x_continuous(limits = c(2034,2100), breaks = seq(2035,2095,10)) +
  scale_size_area(breaks = c(1, 5, 10, 15, 20), limits = c(0, 20)) +  labs(fill = "NTD Change") +
  theme(axis.title = element_text(size=12), legend.position = c("left"), plot.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), #legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8, hjust = 0.5),
      axis.text.y = element_text(hjust = 1, size = 10))
ntdpanel_forpatch


###combo panel####
library(patchwork)
patch <- ntdpanel_forpatch | aapanel_no_yaxislab 
patch<- patch + plot_annotation(title = 'SPBS QS Changes (NTDs on left and AAs on right)', theme = theme(plot.title = element_text(hjust = 0.5)))

pdf("combo_bubbles_spbs_NTDandAA.pdf", width = 14, height = 5)
patch
dev.off()


####heatmap of ABs bubbles for sup fig 14 #####
absite<-read.xlsx("ab_with_srbd.xlsx")

absitelong<-melt(absite, id.vars=c("antibody"))
absitelong
colnames(absitelong)<-c("antibody", "position", "hit")
absitelong$hit
absitelong$hit<-as.character(absitelong$hit)
absitelong$position<-as.numeric(as.character(absitelong$position))
unique(absitelong$antibody)
aborder<-c("P2B-2F6 hmAb", "CB6 hmAb",  "B38 hmAb", "REGN10933", "REGN10987",  "H11-H4 (nAb)", "H11-D4 (nAb)", "C5 (nAb)",  "H3 (nAb)" )
aborder<-rev(aborder)

absitelong$antibody<-factor(absitelong$antibody, levels = aborder)
abcols<-c("0" = "white", "1" = "black")

absitelong_cut<-subset(absitelong, absitelong$hit == 1)

###dots where the dot is scaled to the num of hits
library(dplyr)
absitelong_cut$hitnum<-as.numeric(absitelong_cut$hit)
absitelong_cut_summed <- as.data.frame(absitelong_cut) %>% 
  group_by(position) %>% 
  summarise(hitcount = sum(hitnum))
absitelong_cut_summed$sites<-"all"
absitelong_cut_summed<-subset(absitelong_cut_summed, absitelong_cut_summed$position >471)
absitelong_cut_summed<-subset(absitelong_cut_summed, absitelong_cut_summed$position <505)

abhmdot<-ggplot(absitelong_cut_summed, aes(position, sites, size = hitcount)) + geom_point(alpha=0.75) +  scale_x_continuous(limits = c(470,505), breaks = seq(471,504,1)) +
    theme(axis.title = element_text(size=5), legend.position = c("bottom"), legend.direction = ("horizontal"), plot.background = element_blank(), panel.background = element_blank(),
      plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5), panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      axis.text.x = element_text(size = 6, hjust = 0.5),
      axis.text.y = element_blank())
abhmdot

pdf("absite_map_dots.pdf", width = 7, height = 3)
abhmdot
dev.off()


####end bubbles of ABs for sup fig 14####
