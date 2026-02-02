########################################
### Prepare the data for ulf scripts ###
########################################

#### libraries
library(regioneR)

#### snakemake parameters
# output_file_name = snakemake@output[["done_step"]]
# flair_corrected_bed = snakemake@params[["flair_corrected_bed"]]
# exons_folder = snakemake@parms[["exons_folder"]]

output_file_name = "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/flair/prepare_ulf/sample_done.txt"
flair_corrected_bed = "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/results/flair/correction/LRS4_LRS4_barcode01.flair.corrected_all_corrected.bed"
exons_folder = "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/ref_files/exons/"

#### Do the work
### list all exons df in exons_folder
exons_files <- list.files(path = exons_folder, pattern = "\\.tsv$", full.names = TRUE)
exons_files_list <- lapply(X = exons_files, FUN = function(x) {
  gene_table <- read.table(file = x, header = TRUE, sep = "\t", row.names = NULL)
  colnames(gene_table) <- c("chromosome","start","end","exon_name","exon#","strand") 
  return(gene_table)
})
names(exons_files_list) <- gsub(pattern = "\\.tsv$", replacement = "", x = basename(exons_files))


### Corrected bed file into GRanged object
auxDf <- read.table(file = flair_corrected_bed, header = FALSE, stringsAsFactors = FALSE)
correctedRanges = regioneR::toGRanges(auxDf[,c(1:3)])

for (i in 1:length(names(exons_files_list))) {
  ## gene now 
  gene_now <- names(exons_files_list)[i]
  exons_now <- exons_files_list[[i]]
  
  ## final_name
  sample_name <- gsub(pattern = "\\.flair\\.corrected_all_corrected\\.bed$", 
                      replacement = "", x = basename(flair_corrected_bed))
  final_place <- gsub(pattern = "sample_done\\.txt$", 
                      replacement = "", x = output_file_name)
  final_file_name <- paste0(final_place,sample_name,gene_now,".all_corrected.overlapping_BRCA1.bed")
  findover
}

# auxDF <- read.table(flair_corrected_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# correctedRanges = regioneR::toGRanges(auxDF[,c(1:3)])

#### write the output
cmd <- paste0("touch ",output_file_name)
system(cmd)