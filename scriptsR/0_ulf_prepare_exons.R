#################################################
### Filter transcripts out of a transcriptome ###
#################################################

#### libraries
library(biomaRt)

#### snakemake parameters
genes_to_filter <- snakemake@params[["genes_to_filter"]]
place_of_final_information <- snakemake@params[["place_of_final_information"]]
### for now
# genes_to_filter <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this_ENIGMA.txt"
# place_of_final_information <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/ref_files/exons/"

genes_to_filter <- readLines(genes_to_filter)

#### Connect to biomaRt and do the job
mart <- biomaRt::useEnsembl(biomart = "ensembl", version = 115)
datasets <- biomaRt::listDatasets(mart = mart)
ensembl <- biomaRt::useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")

### get the filters
filts <- biomaRt::listFilters(mart = ensembl)
attrs <- biomaRt::listAttributes(mart = ensembl)

annotation_of_genes_select <- biomaRt::getBM(mart = ensembl,
                                             attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id","transcript_mane_select"),
                                             filters = c("hgnc_symbol"), 
                                             values = genes_to_filter)

annotation_of_genes_select <- annotation_of_genes_select[annotation_of_genes_select[[4]] != "",, 
                                                         drop = FALSE]
# annotation_of_genes_select <- annotation_of_genes_select[[2]]

annotation_of_genes <- biomaRt::getBM(mart = ensembl,
                                      attributes = c("ensembl_transcript_id",
                                                     "chromosome_name",
                                                     "exon_chrom_start","exon_chrom_end",
                                                     "ensembl_exon_id",
                                                     "rank","strand"),
                                      filters = c("ensembl_transcript_id"), 
                                      values = annotation_of_genes_select[[3]])

### generate annotation_final
annotation_final <- merge(annotation_of_genes_select, annotation_of_genes, by = "ensembl_transcript_id")
cols_i_want <- c("hgnc_symbol","chromosome_name","exon_chrom_start","exon_chrom_end",
                 "ensembl_exon_id","rank","strand")
annotation_final <- annotation_final[,grep(pattern = paste0(x = cols_i_want, collapse = "|"), x = colnames(annotation_final)), drop = FALSE]
final_col_names <- c("hgnc_symbol","chromosome","start","end","exon_name","exon#","strand")
colnames(annotation_final) <- final_col_names

### save the final annotations per gene
genes_to_filter[which(!genes_to_filter %in% unique(annotation_final$hgnc_symbol))]
genes_to_filter <- genes_to_filter[genes_to_filter %in% unique(annotation_final$hgnc_symbol)]
for (i in 1:length(genes_to_filter)) {
  ## data now
  gene_now <- genes_to_filter[i]
  annotation_now <- annotation_final[annotation_final[[1]] == gene_now,, drop = FALSE]
  ## strand name change
  strand_old <-unique(annotation_now$strand)
  strand_new <- ""
  if (strand_old == -1) {
    strand_new <- "-"
  } else if (strand_old == 1) {
    strand_new <- "+"
    }
  annotation_now$strand <- strand_new
  ## exon rank order
  annotation_now <- annotation_now[order(annotation_now[[6]]),, drop = FALSE]
  ## gene removal
  annotation_now <- annotation_now[,-1, drop = FALSE]
  ## save the file
  write.table(x = annotation_now, file = paste0(place_of_final_information,gene_now,".tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  }


### file all done
cmd <- paste0("touch ",place_of_final_information,"AAA_done.txt")
system(cmd)
