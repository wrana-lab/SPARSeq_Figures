R scripts for Figure 6 and examples of Supplemental Figure 12 analysis:

After developing our quasispecies analysis pipeline available at https://github.com/wrana-lab/SPARSEQ_QUASISPECIES we needed to develop reasomable filtering methods to determine putative quasispecies sequences (pQS) vs. background noise sequences. This folder includes an example of the processing method using the Delta samples. The same method was followed for all subsets of samples.

Figure 6 A shows the flowchart of our analysis process.

Figure 6 B shows the finalized count of S-Rbm pQS.

Figure 6 C shows the finalized aligned S-Rbm pQS per VOC.

Supplemental Figure 12 shows the development of the pQS analysis as shown in the flowchart in Fig. 6 A: 

Sup. Fig. 12 A shows the distributions of S-Rbm counts across samples per barcoded run.

Sup. Fig. 12 B shows the series of test cutoffs and linear models we used for our S-Rbm data, to determine a finalized, conservative cutoff for our quasispecies analysis for each run.

Sup. Fig. 12 C shows the same analysis as in B, but for the S-Pbs data. 

Sup. Fig. 12 D shows the finalized count of S-Pbs pQS.

Sup Fig 12 E shows the aligned S-Pbs pQS per VOC.


The script ExampleSubsetProcessing.R shows an example of how we compiled the Delta subset quasispecies runs at different percentage cutoffs to determine a finalized cutoff. 

The script ExampleFinalizedCutoffProcessing.R shows an example of how we compare the two final copies of the bc1 and bc2 data for our Delta subset, using the x-intercept based cutoffs from the above example, to get a final list of pQS which are seen in both copies of individual samples. In this script the Bioconductor package Biostrings (https://bioconductor.org/packages/release/bioc/html/Biostrings.html) is needed for easy translation of nucleotides to amino acids.

The outputs of the second script are later aggregated per variant - see subfolder Figure6CandSupplemental13-14.

For SPbs analysis, we did the same workflow but with a modified version of the quasispecies pipeline that didn't use barcodes, since we didn't use barcodes on the SPbs amplicons. 
