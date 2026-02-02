#################################################
### Filter transcripts out of a transcriptome ###
#################################################

#### libraries
library(seqinr)

#### The genes of interest
# genes_of_int <- c("LZTR1","NF2","PKD2","SMARCE1","TSC1","NF1","PKD1","SMARCB1","SPRED1","TSC2")
genes_of_int = snakemake@params[["genes_to_work_on"]]

#### get the transcriptome
# transcriptome <- "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/gencode.v48.pc_transcripts.fa"
transcriptome = snakemake@input[["transcriptome"]]
transcriptome = read.fasta(transcriptome, seqtype = "DNA", forceDNAtolower = FALSE)
out_file = snakemake@output[["transcriptome_filtered"]]

for (i in 1:length(genes_of_int)) {
  trans_to_keep <- names(transcriptome)[grepl(pattern = paste0("\\|",genes_of_int[i],"\\|"), x = names(transcriptome))]
  transcriptome_now <- transcriptome[trans_to_keep]
  write.fasta(
    transcriptome_now, names = trans_to_keep,
    file.out = paste0(out_file), as.string = FALSE)        
#   write.fasta(transcriptome_now, names = trans_to_keep,   
#   file.out = paste0("/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/per_gene/",genes_of_int[i],".fasta"), as.string = FALSE)
}

