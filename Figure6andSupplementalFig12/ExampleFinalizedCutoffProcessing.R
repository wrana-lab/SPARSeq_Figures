##This script uses our Delta data as an example to show how we process the final version of our quasispecies pipeline,
#using the x-intercepts we calculated from the previous script, to get a final list of pQS to move forward with analysis on.



###Set up reference sequences so we can easily compare later in a visualization tool 
srbd_omicron<-"ATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT"
srbd_omicron_AA<-"IYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGV"
srbd_delta<-"ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"
srbd_delta_AA<-"IYQAGSKPCNGVEGFNCYFPLQSYGFQPTNGV"
srbd_wt<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"
srbd_wt_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGV"
srbd_alpha<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTTATGGTGTT"
srbd_alpha_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGV"


#BiocManager::install("Biostrings")
library(Biostrings)
### Translating DNA/RNA:
#translate(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

#import aligned lists from the final, x-intercept based analysis runs of the quasispecies pipeline: Cutoff of copy 1: 0.2511, Cutoff of copy 2: 0.2027
bc1_aligned_list<-read.csv("bc1_2022/srbd_aligned_list.txt")
bc2_aligned_list<-read.csv("bc2_2023/srbd_aligned_list.txt")

#drop refseq from first row 
bc1_aligned_list<-bc1_aligned_list[2:nrow(bc1_aligned_list),]
bc2_aligned_list<-bc2_aligned_list[2:nrow(bc2_aligned_list),]

nrow(bc1_aligned_list) #583
nrow(bc2_aligned_list) #322

#remove gaps from the aligned sequences etc 
bc1_aligned_list$sampleID<-substr(bc1_aligned_list$sample, 1, 9)
bc1_aligned_list$seq_chars<-gsub("-", "", bc1_aligned_list$sequence)
bc1_aligned_list$seq_chars<-gsub("C$", "", bc1_aligned_list$seq_chars)

bc2_aligned_list$sampleID<-substr(bc2_aligned_list$sample, 1, 9)
bc2_aligned_list$seq_chars<-gsub("-", "", bc2_aligned_list$sequence)
bc2_aligned_list$seq_chars<-gsub("C$", "", bc2_aligned_list$seq_chars)

unique(bc1_aligned_list$sampleID) #112 kept samples
unique(bc2_aligned_list$sampleID) # 135

bc1_aligned_list<-bc1_aligned_list[,3:4]
bc2_aligned_list<-bc2_aligned_list[,3:4]

#check count output list - this is the output from the quasispecies pipeline that just contains total Srbm counts per sample
bc1countlist<-read.csv("bc1_2022/countoutputlist.csv")
bc2countlist<-read.csv("bc2_2023/countoutputlist.csv")

bc1countlist<-bc1countlist[bc1countlist$total_SRBD_count > 32000,]
bc2countlist<-bc2countlist[bc2countlist$total_SRBD_count > 32000,]
#keep only those samples which had at least 32k detected Srbm reads according to the python scipt analysis.

##combine the sample ID and sequence because we need to be able to compare if the same sample has the same sequences in both copies of the analysis 
bc1_aligned_list_merged<-bc1_aligned_list #bc1_aligned_list_merged$exp<-NULL
bc2_aligned_list_merged<-bc2_aligned_list #bc2_aligned_list_merged$exp<-NULL
bc1_aligned_list_merged$mergedsamplesequence<-paste(bc1_aligned_list_merged$sampleID, bc1_aligned_list_merged$seq_chars, sep = "_")
bc2_aligned_list_merged$mergedsamplesequence<-paste(bc2_aligned_list_merged$sampleID, bc2_aligned_list_merged$seq_chars, sep = "_")

#we only want to keep those sequences which are found in both copies of a sample
bc1_final<-bc1_aligned_list_merged[bc1_aligned_list_merged$mergedsamplesequence %in% bc2_aligned_list_merged$mergedsamplesequence,]
bc2_final<-bc2_aligned_list_merged[bc2_aligned_list_merged$mergedsamplesequence %in% bc1_aligned_list_merged$mergedsamplesequence,]


#create function  table of common sample / seq
compile_common_qs <- function(x,y){
  common_aligned_list_sm <- x[0,]
  for(i in unique(y$sampleID)){
    sub_bc1<-subset(x, x$sampleID == i)
    sub_bc2<-subset(y, y$sampleID == i)
    for(j in unique(sub_bc2$seq_chars)){
      if (any(grepl(j, sub_bc1$seq_chars, fixed=TRUE))){
        #push sample id and sequence to output table 
        common_aligned_list_sm[nrow(common_aligned_list_sm)+1,] <- c(i, j)
      }
    }
  }
  return(common_aligned_list_sm)
}

bc1_aligned_list_sm<-bc1_final[,1:2]
bc2_aligned_list_sm<-bc2_final[,1:2]
bcmin1<-compile_common_qs(bc1_aligned_list_sm, bc2_aligned_list_sm)

bcmin1_agg<-bcmin1 %>% 
  group_by(seq_chars) %>% 
  summarise(num_samples = n())
bcmin1_agg<-bcmin1_agg[order(bcmin1_agg$num_samples, decreasing = T),]
bcmin1_agg<-as.data.frame(bcmin1_agg)
nrow(bcmin1_agg) #
common_aligned_list_sm<-bcmin1_agg

write.xlsx(bcmin1, file = "delta_pQS_sample_sequence_list.xlsx")
unique(bcmin1$sampleID) #12

write.xlsx(common_aligned_list_sm, "delta_pQS_xint_cutoffs_common_sampleIDsandsequences_finalized.xlsx")

#convert to AAs
common_dnastring<-DNAStringSet(common_aligned_list_sm$seq_chars)
common_aas<-as.data.frame(translate(common_dnastring))
common_translated<-cbind(common_aas,common_aligned_list_sm)
colnames(common_translated)<-c("aas", "sequence", "numSamples")


###next we write out the sequences in fasta format with multiple references sequences
##we do this so we can easily open these in visualization software and easily see the pQS sequences with
#context of all main VOC sequences 
write.table(">Refseq_srbd_omicron", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = FALSE, col.names = FALSE, quote = F)
write.table(srbd_omicron, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_delta", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_delta, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_WT", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_wt, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_alpha", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_alpha, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)

for(i in 1:nrow(common_translated)){
  seq<-common_translated[i,2]
  exp_samp<-common_translated[i,3]
  lines<-paste0(">",exp_samp,"_samples")
  write.table(lines, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
  write.table(seq, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
}

#AAs
write.table(">Refseq_srbd_omicron", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = FALSE, col.names = FALSE, quote = F)
write.table(srbd_omicron_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_delta", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_delta_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_WT", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_wt_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_alpha", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_wt_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)

for(i in 1:nrow(common_translated)){
  seq<-common_translated[i,1]
  exp_samp<-common_translated[i,3]
  lines<-paste0(">",exp_samp,"_samples")
  write.table(lines, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
  write.table(seq, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
}

#To get an aligned fasta file, we use clustal omega installed locally and run:
# clustalo -i delta_srbd_topAAsequences.fa -o delta_srbd_topAAsequences.fasta
#for each of the results files.


