#######################
### GFF3 filtration ###
#######################

### libraries
library(readr)
library(dplyr)

### Snakemake variables
raw_gff3 <- snakemake@input[["gff3_filtered"]]
final_gff3 <- snakemake@output[["gff3_filtered_MANE"]]
# raw_gff3 <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/ref_files/filtered.gff3"
# final_gff3 <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/ref_files/filtered_MANE2.gff3"

gff <- readr::read_tsv(file = raw_gff3, skip = 7, col_names = FALSE)
gff3_header <- readLines(con = raw_gff3, n = 7)

colnames(gff) <- c("chr","source_of_ann","type","start","end","dot1","strand","dunno","tag")
valid_chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                "chr21","chr22","chrX","chrY") 

### filter he gff3 file
gff_filtered <- gff %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(keep = ifelse(type == "gene","YES, gene",
                              ifelse(grepl(x = tag, pattern = "Mane_Select|MANE_Select"),"YES, MANE","NO"))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(keep %in% c("YES, MANE","YES, gene")) %>% 
  dplyr::filter(chr %in% valid_chrs) %>% 
  dplyr::select(-c(keep))

### check the correct filtration of the gene's data
## which genes are being checked
genes_to_check <- gff_filtered$tag[gff_filtered$type == "gene"]
genes_to_check <- gsub(pattern = "\\..*", replacement = "", x = genes_to_check)
genes_to_check <- gsub(pattern = "ID\\=", replacement = "", x = genes_to_check)

## get the gene - transcripts from the gff3 files
gff_check <- gff_filtered %>%  
  dplyr::filter(type %in% c("gene","transcript"))

bad_genes <- c()
for (i in 1:length(genes_to_check)) {
  gff_check_now <- gff_check %>% 
    dplyr::filter(grepl(x = tag, pattern = genes_to_check[i]))
  if (nrow(gff_check_now) == 2) {
    next
  } else if (nrow(gff_check_now) != 2) {
    bad_genes <- c(bad_genes,genes_to_check[i])
  }
  
}

gff_filtered <- gff_filtered %>%
  dplyr::filter(!grepl(pattern = paste0(x = c(bad_genes, collapse = "|")), x = tag))

bag_genes_df <- data.frame()
if (length(bad_genes) != 0) {
  for (i in 1:length(bad_genes)) {
    ## Gene information
    bad_genes_now <- bad_genes[i]
    gff_gene <- gff %>% 
      dplyr::filter(type == "gene") %>% 
      dplyr::filter(grepl(x = tag, pattern = bad_genes_now))
    
    ## Transcript information
    gff_transcripts <- gff %>% 
      dplyr::filter(type == "transcript") %>% 
      dplyr::filter(grepl(x = tag, pattern = bad_genes_now)) %>% 
      dplyr::mutate(length_of_t = abs(start - end))
    gff_transcripts <- gff_transcripts %>% 
      dplyr::filter(length_of_t == max(length_of_t, na.rm = TRUE)) %>% 
      dplyr::select(-length_of_t)
    gff_transcripts <- gff_transcripts[1,, drop = FALSE]
    gff_transcripts_id <- gff_transcripts$tag 
    gff_transcripts_id <- gsub(pattern = "\\..*", replacement = "", x = gff_transcripts_id)
    gff_transcripts_id <- gsub(pattern = "ID\\=", replacement = "", x = gff_transcripts_id)
    
    ## All other information
    other_info <- gff %>% 
      dplyr::filter(!type %in% c("gene","transcript")) %>% 
      dplyr::filter(grepl(x = tag, pattern = paste0(x = c("Parent=",gff_transcripts_id), collapse = "")))
    
    list_to_bind <- list(gff_gene,
                         gff_transcripts,
                         other_info)
    
    bad_genes_info <- do.call(rbind, list_to_bind)
    bad_genes_info <- as.data.frame(bad_genes_info)
    bag_genes_df <- rbind(bag_genes_df,bad_genes_info)
    }
  gff_filtered <- rbind(gff_filtered,bad_genes_info)
  }

writeLines(gff3_header, con = final_gff3)
readr::write_tsv(x = gff_filtered, 
                 file = final_gff3, 
                 append = TRUE, col_names = FALSE)

readr::write_tsv(file = final_gff3,
                 x = gff_filtered, 
                 col_names = FALSE)
