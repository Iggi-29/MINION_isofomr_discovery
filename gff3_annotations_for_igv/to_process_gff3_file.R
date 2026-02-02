#####################################################
### Get gff3 file and get MANE Select Transcripts ###
#####################################################

### libraries
library(tidyverse)
library(dplyr)
library(readr)

### WD
setwd("/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION")

### gff3 file
gff3_annotation_file <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gff3"

## body
gff3_annotation <- readr::read_tsv(
  file = gff3_annotation_file,
  skip = 7, col_names = FALSE)
colnames(gff3_annotation) <- c(
  "chr","source_of_ann",
  "type","start","end",
  "dot1","strand","dunno","tag")

## header
gff_header <- readLines(con = gff3_annotation_file, n = 7)

# gff3_annotation <- gff3_annotation %>% 
#   dplyr::filter(chr == "chr22")

### Filter for genes and MANE Select / MANE Select Clinical objects
gff3_annotation <- gff3_annotation %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(keep = ifelse(type == "gene","YES, gene","NO")) %>% 
  dplyr::mutate(keep = ifelse(grepl(x = tag, pattern = "Mane|MANE"),"YES, Mane",keep)) %>% 
  dplyr::ungroup()

gff3_annotation <- gff3_annotation %>% 
  dplyr::filter(grepl(x = keep, pattern = "YES")) %>% 
  dplyr::select(-c(keep))


### Change ID in ensemble to the gene name to see it in IGV 
for (i in 1:length(rownames(gff3_annotation))) {
  ## data now
  if (gff3_annotation$type[i] == "gene") {
    tag_now <- gff3_annotation$tag[i]
    id <- gsub(pattern = ".*ID=([^;]+);.*", replacement = "\\1", x = tag_now)
    gene <- gsub(pattern = ".*gene_name=([^;]+);.*", replacement = "\\1", x = tag_now)
    tag_now <- gsub(pattern = "ID=[^;]+", replacement = paste0("ID=", gene), x = tag_now)
    
    gff3_annotation$tag[i] <- tag_now
    cat(paste0("Done with row: ",i," out of ",length(rownames(gff3_annotation)),"\r"))} else if (gff3_annotation$type[i] != "gene") {
      next
      cat(paste0("Done with row: ",i," out of ",length(rownames(gff3_annotation)),"\r"))}
  }

### write the file
# writeLines(
#   text = gff_header,
#   con = "./gff3_annotations_for_igv/gff3_filtered.gff3")
# 
# readr::write_tsv(
#   file = "./gff3_annotations_for_igv/gff3_filtered.gff3",
#   x = gff3_annotation, col_names = FALSE, append = TRUE)

