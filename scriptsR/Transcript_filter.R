#################################################
### Filter transcripts out of a transcriptome ###
#################################################

#### libraries
library(seqinr)

#### get the transcriptome
transcriptome <- "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/gencode.v48.pc_transcripts.fa"
transcriptome <- read.fasta(transcriptome, seqtype = "DNA", forceDNAtolower = FALSE)

genes_of_int <- c("LZTR1","NF2","PKD2","SMARCE1","TSC1","NF1","PKD1","SMARCB1","SPRED1","TSC2")
for (i in 1:length(genes_of_int)) {
  trans_to_keep <- names(transcriptome)[grepl(pattern = paste0("\\|",genes_of_int[i],"\\|"), x = names(transcriptome))]
  transcriptome_now <- transcriptome[trans_to_keep]
  write.fasta(transcriptome_now, names = trans_to_keep, file.out = paste0("/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/per_gene/",
                                                                          genes_of_int[i],".fasta"), as.string = FALSE)
}

