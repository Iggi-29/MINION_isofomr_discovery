#######################################
### FLAMES long reads data analysis ###
#######################################

#### libraries
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

#### Constants
### save the data here
place_to_save_the_data <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/"
### base path for flame results
flame_base_path_results <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/" ## here they are the same but it could be different
### annotation of the genes
raw_annotation_gft_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/gff3_annotations_for_igv/gff3_filtered.gff3" ## MANE!
raw_annotation_gft_file <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/ref_files/filtered.gff3"
### genes to work on
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this.txt"
### matching final df
matching_final_df <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matching_to_annotate.txt"
### raw annotation gtf
annotation_of_trans_gtf <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/annotation_of_transcripts.xlsx"
### output file exon
output_file_exon <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matched_exon_data"
### output file genes and trans
output_file_genes_and_trans <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matched_gene_and_trans_data"

######### HERE PUT THE SNAKEMAKE PARAMS SO ALL THAT IS ABOVE THIS PART SHOULD BE TAKEN BY SNAKEMAKE
######
###

##### Prepare the script to work on
#### places to save the data
if(!dir.exists(paths = paste0(place_to_save_the_data,"/final"))) {
  dir.create(path = paste0(place_to_save_the_data,"/final"))
  cat(paste0("Created the final data dir","\n"))
} else {
  cat(paste0("Didn't need to create the final data dir","\n"))
}

#### flame results importation
flame_results <- list.files(path = flame_base_path_results,
                            recursive = TRUE, full.names = TRUE,
                            pattern = "isoform_annotated\\.filtered\\.gff3$")
# flame_results <- flame_results[grep(x = flame_results, pattern = "LRS4_LRS4_barcode01|LRS4_LRS4_barcode02")] ## només amb 2 arxius com a prova

#### Genes to work-on
genes_to_work_on <- readLines(con = genes_to_work_on)

#### GTF-like annotation
gtf_like_annotation <- openxlsx2::read_xlsx(file = annotation_of_trans_gtf)

gene_annotation <- gtf_like_annotation %>% 
  dplyr::select(gene_name, gene_id, transcript_id) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(reff_transcript = transcript_id)

gtf_like_annotation <- gtf_like_annotation %>% 
  dplyr::select(c(start,end,ID)) %>% 
  dplyr::rename(original_start = start,
                original_end = end,
                raw_id = ID) %>% 
  dplyr::filter(!duplicated(raw_id)) 

### Get all the flames results too
flame_results_list <- list()
for (i in 1:length(flame_results)) {
  ### Data now
    flame_results_now <- flame_results[i]
    flame_results_now <- rtracklayer::import(con = flame_results_now)
    flame_results_now <- as.data.frame(flame_results_now)
    flame_results_now$Parent <- sapply(
      X = flame_results_now$Parent,
      FUN = function(x) {
        parent_worked_col <- ""
        if (length(x) == 0) {
          parent_worked_col <- NA_character_
        } else if (length(x) != 0) {
          parent_worked_col <- x[1]
        }
        return(parent_worked_col)
        })
    sample_name <- gsub(
      pattern = ".*/(.*)/[^/]+$",
      replacement = "\\1",
      x = flame_results[i])
    flame_results_list[[sample_name]] <- flame_results_now
    }

remove(flame_results_now)

matching_data_list <- list()
flame_gene_and_transcript_data <- list()

#### Matthing df work
matching_final_df <- read.table(file = matching_final_df, header = TRUE, sep = "\t")
for (i in 1:length(unique(matching_final_df$sample_name))) {
  ### Data now
  sample_now <- unique(matching_final_df$sample_name)[i]
  matching_data_now <- matching_final_df %>% 
    dplyr::filter(sample_name == sample_now)
  flame_results_now <- flame_results_list[[sample_now]]
  
  ### Filter for the transcript and gene apps in the flame_results
  flame_results_now_gene_transcript <- flame_results_now %>% 
    dplyr::filter(type %in% c("gene","transcript")) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(flame_type = type)
  
  flame_results_now_exons <- flame_results_now %>% 
    dplyr::filter(!type %in% c("gene","transcript")) %>% 
    dplyr::rename(flame_type = type,
                  flame_transcript = Parent)
  
  ### Get the info inside the matching data now
  matching_data_now <- matching_data_now %>%
    dplyr::left_join(y = flame_results_now_exons, by = c("exon_id","flame_transcript")) %>% 
    dplyr::select(-c(score,phase,ID,gene_id,transcript_id)) %>% 
    dplyr::rename(transcript_id = flame_transcript) %>% 
    dplyr::mutate(Parent = transcript_id) %>% 
    dplyr::mutate(flame_type = "exon")
  
  ### Get the gene annotation required data
  matching_data_now <- matching_data_now %>% 
    dplyr::left_join(y = gtf_like_annotation, by = "raw_id") %>% 
    dplyr::mutate(diff_start = abs(start - original_start),
                  diff_end = abs(end - original_end)) %>% 
    dplyr::left_join(y = gene_annotation, by = "reff_transcript")
  
  ## Exon range
  matching_data_now <- matching_data_now %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::mutate(exon_range = {
      max_range <- max(as.numeric(rank), na.rm = TRUE)
      min_range <- min(as.numeric(rank), na.rm = TRUE)
      
      paste0(min_range,"-",max_range)
    })

  ### id
  matching_data_now <- matching_data_now %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(id = paste0(sample_name,"_",gene_name,"_",transcript_id)) %>% 
    dplyr::ungroup()
  
  ### raw_id
  matching_data_now <- matching_data_now %>% 
    dplyr::rename(ensembl_coincidence_id = raw_id) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(raw_id = {
      before <- gsub(pattern = ":.*$", replacement = "", x = ensembl_coincidence_id)
      after <- gsub(pattern = "^.*:", replacement = "", x = ensembl_coincidence_id)
      
      ifelse(type %in% c("exon","intron"),
             paste0(before,"_",after),type)
      
    }) %>% 
    dplyr::ungroup()
  
  ### finish the dataframe
  matching_data_now <- matching_data_now %>% 
    dplyr::rename(ensembl_gene_id = gene_id,
                  gene_symbol = gene_name,
                  info_from = flame_type) %>% 
    dplyr::select(sample_name, seqnames, start, end, original_start, original_end, diff_start, diff_end,
                  strand, width, source, type, exon_id, gene_symbol, ensembl_gene_id, 
                  transcript_id, reff_transcript,raw_id ,ensembl_coincidence_id, exon_range, Parent, 
                  rank, overlap_width, percent_overlap, info_from,
                  support_count, id)
  
  ### match the gene and trans data to the exon data
  flame_results_now_gene_transcript <- flame_results_now_gene_transcript %>% 
    dplyr::select(-c(score, phase)) %>% 
    dplyr::rename(raw_id = ID,
                  ensembl_gene_id = gene_id,
                  info_from = flame_type) %>% 
    dplyr::mutate(sample_name = sample_now)
  # flame_results_now_gene_transcript <- flame_results_now_gene_transcript %>% 
  #   dplyr::
  ### put the results into a tabular format
  matching_data_list[[sample_now]] <- matching_data_now
  flame_gene_and_transcript_data[[sample_now]] <- flame_results_now_gene_transcript
  }

#### Write the output
### Exon matched data
openxlsx::write.xlsx(x = matching_data_list,
                     file = paste0(output_file_exon,".xlsx"))
readr::write_tsv(x = do.call(rbind ,matching_data_list),
                     file = paste0(output_file_exon,".tsv"))
### Gene and transcript information
openxlsx::write.xlsx(x = flame_gene_and_transcript_data,
                     file = paste0(output_file_genes_and_trans,".xlsx"))
readr::write_tsv(x = do.call(rbind,flame_gene_and_transcript_data),
                     file = paste0(output_file_genes_and_trans,".tsv"))




