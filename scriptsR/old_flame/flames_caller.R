#######################################
### FLAMES long reads data analysis ###
#######################################

#### libraries
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(ggnewscale)
library(ggrepel)

### functions
## ranges generator
ranges_generator <- function(gff_transformed_df) {
  ranges_list <- list()
  for (i in 1:length(unique(gff_transformed_df$transcript_id))) {
    ## data now
    transcript_now <- unique(gff_transformed_df$transcript_id)[i]
    gene_now <- unique(gff_transformed_df$gene_name[gff_transformed_df$transcript_id == transcript_now])
    df_now <- gff_transformed_df %>% 
      dplyr::filter(transcript_id == transcript_now) %>% 
      dplyr::filter(!type %in% "transcript")
    
    df_now_cds <- df_now %>% 
      dplyr::filter(type == "CDS") %>% 
      dplyr::filter(exon_number == max(as.numeric(exon_number), na.rm = TRUE) |
                      exon_number == min(as.numeric(exon_number), na.rm = TRUE))
    
    df_now <- df_now %>% 
      dplyr::filter(!type %in% c("CDS"))
    
    df_now <- rbind(df_now, df_now_cds)
    
    ## raw_id generation
    df_now <- df_now %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(raw_id = paste0(as.character(type),"_",as.character(exon_number))) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(start, end, raw_id)
    
    ## id of the ranges
    range_id <- paste0(gene_now,"_",transcript_now)
    
    ## Generate IRanges data
    ranges_now <- IRanges::IRanges(start = df_now$start, 
                                   end = df_now$end, 
                                   names = df_now$raw_id)
    
    ranges_list[[range_id]] <- ranges_now
    
  }
  
  return(ranges_list)
}
## comparison maker
comparison_maker <- function(flame_results_df, ranges_to_check) {
  df_comparisons <- data.frame()
  for (i in 1:length(unique(flame_results_df$id))) {
    ### data now
    id_now <- unique(flame_results_df$id)[i]
    data_now <- flame_results_df %>% 
      dplyr::filter(id == id_now)
    ## info to annotate
    gene_now <- data_now$gene_symbol[data_now$type == "gene"]
    trans_now <- data_now$transcript[data_now$type == "transcript"]
    trans_now <- gsub(pattern = "transcript\\:", replacement = "", x = trans_now)
    sample_name_now <- unique(data_now$sample_name)
    ## def data_now
    data_now <- data_now %>%
      dplyr::filter(!type %in% c("gene","transcript"))
    
    ## ranges now
    ranges_gene_now <- ranges_to_check[grep(pattern = gene_now, x = names(ranges_to_check))]
    
    ### from data_now generate IRanges
    ranges_now <- IRanges::IRanges(start = data_now$start,
                                   end = data_now$end, 
                                   names = data_now$raw_id)
    
    ### Do the intersections
    for (app in 1:length(ranges_gene_now)) {
      name_of_query <- names(ranges_gene_now)[app]
      ranges_gene_now_trans <- ranges_gene_now[[app]]
      
      hits <- findOverlaps(subject = ranges_now,
                           query = ranges_gene_now_trans,
                           minoverlap = 0)
      
      overlaps <- pintersect(y = ranges_now[subjectHits(hits)],
                             x = ranges_gene_now_trans[queryHits(hits)])
      
      percent_overlap <- width(overlaps)/width(ranges_gene_now_trans[queryHits(hits)])
      
      matching_df <- data.frame(
        subject_id = names(ranges_now[subjectHits(hits)]),
        refference_id = names(ranges_gene_now_trans[queryHits(hits)]),
        overlap_width = width(overlaps),
        percent_overlap = round(x = percent_overlap, digits = 4),
        name_of_query = name_of_query)
      
      ### Finish matching df
      matching_df <- matching_df %>% 
        dplyr::rowwise() %>% 
        dplyr::filter(subject_id == refference_id) %>% 
        dplyr::ungroup()
      
      ### check for 0 overlapign parts of the transcript
      all_parts_of_trans <- names(ranges_gene_now_trans)
      no_present_all_parts_of_trans <- names(ranges_gene_now_trans)[!all_parts_of_trans %in% matching_df$refference_id]
      
      if (length(no_present_all_parts_of_trans) == 0) {
        matching_df <- matching_df
      } else if (length(no_present_all_parts_of_trans) != 0) {
        df_to_bind <- data.frame(
          subject_id = "NO present",
          refference_id = no_present_all_parts_of_trans,
          overlap_width = 0,
          percent_overlap = 0,
          name_of_query = name_of_query)
        
        matching_df <- rbind(matching_df, df_to_bind)
      }
      matching_df$gene <- gene_now
      matching_df$transcript <- trans_now
      matching_df$sample_name <- sample_name_now
      matching_df$id <- id_now
      
      ### finish the data.frames
      df_comparisons <- rbind(df_comparisons,matching_df)
      
    }
  }
  return(df_comparisons)
}
## add introns to a gtf_df
add_introns_to_gtf <- function(gtf_to_add_info) { 
  ## save the original data 
  gtf_to_add_info <- gtf_to_add_info 
  ## calculate the introns 
  gtf_intron <- gtf_to_add_info %>% 
    ## do it with introns 
    dplyr::filter(type == "exon") %>% 
    ## by each transcript 
    dplyr::group_by(transcript_id) %>% 
    ## do the calculation 
    dplyr::arrange(start, .by_group = TRUE) %>% 
    dplyr::mutate(intron_start = lag(end+1), 
                  intron_end = start -1) %>% 
    dplyr::ungroup() %>% 
    ## remove non valid introns 
    dplyr::filter(!is.na(intron_start)) %>% 
    ## change colnames for a clean rbind 
    dplyr::relocate(intron_start, .after = start) %>% 
    dplyr::relocate(intron_end, .after = end) %>% 
    dplyr::select(-c(start, end)) %>% 
    dplyr::rename(start = intron_start, end = intron_end) %>% 
    dplyr::mutate(type = "intron") 
  gtf_intron <<- gtf_intron
  
  gtf_to_add_info <- rbind(gtf_to_add_info, gtf_intron) 
  return(gtf_to_add_info)
  }
## filter CDS
filter_CDS <- function(raw_annotation, exon_n_col, groupping) {
  raw_annotation <- raw_annotation
  ## filter max and min CDS
  raw_annotation_cds <- raw_annotation %>% 
    dplyr::filter(type == "CDS") %>% 
    dplyr::group_by(.data[[groupping]]) %>% 
    dplyr::filter(.data[[exon_n_col]] == max(as.numeric(.data[[exon_n_col]]), na.rm = TRUE) | 
                    .data[[exon_n_col]] == min(as.numeric(.data[[exon_n_col]]), na.rm = TRUE)) %>% 
    dplyr::ungroup()
  ## filter non CDS
  raw_annotation_temp <- raw_annotation %>% 
    dplyr::filter(type != "CDS")
  ## do the final CDS
  final_annot <- rbind(raw_annotation_temp, raw_annotation_cds)
  return(final_annot)
}

#######################
### Raw data import ###
#######################

#### flame results
flame_results_file_name <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flames/final/annotated_flames_data.xlsx"
flame_results <- openxlsx::read.xlsx(xlsxFile = flame_results_file_name)
flame_results <- flame_results %>% 
  dplyr::filter(sample_name == "LRS4_LRS4_barcode01") %>% 
  dplyr::filter(gene_symbol %in% c("MSH6","PMS2","APC"))  %>% 
  dplyr::filter(grepl(x = transcript, pattern = "ENST"))
# flame_results <- flame_results %>% 
#   # dplyr::filter(transcript == "transcript:ENSG00000134982.19_112737885_112791833_1") %>% 
#   dplyr::filter(sample_name == "LRS4_LRS4_barcode01") %>% 
#   dplyr::filter(gene_symbol %in% c("BRCA1","APC"))

#### genes to work on
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this_ENIGMA.txt"

#### Biomart data importation
genes_to_work_on <- readLines(con = genes_to_work_on)
genes_to_work_on <- c(genes_to_work_on)
### connect to the service
mart <- biomaRt::useEnsembl(biomart = "ensembl", version = 115)
datasets <- biomaRt::listDatasets(mart = mart)
ensembl <- biomaRt::useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
### get the filters
filts <- biomaRt::listFilters(mart = ensembl)
attrs <- biomaRt::listAttributes(mart = ensembl)
### get the annotations
annotation_of_genes <- biomaRt::getBM(mart = ensembl,
                                      attributes = c("ensembl_gene_id","ensembl_gene_id_version",
                                                     "hgnc_symbol","ensembl_transcript_id_version","ensembl_exon_id"),
                                      filters = c("hgnc_symbol"),
                                      values = genes_to_work_on)
annotation_of_genes_select <- biomaRt::getBM(mart = ensembl,
                                             attributes = c("hgnc_symbol","ensembl_transcript_id_version","transcript_mane_select"),
                                             filters = c("hgnc_symbol"), values = genes_to_work_on)
annotation_of_genes_select <- annotation_of_genes_select %>%
  dplyr::filter(transcript_mane_select != "")  

#### GTF annotation Select
raw_annotation_gft_file_select <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/gff3_annotations_for_igv/gff3_filtered.gff3"
raw_annotation_gtf_select <- rtracklayer::import(con = raw_annotation_gft_file_select)
raw_annotation_gtf_select <- raw_annotation_gtf_select[mcols(raw_annotation_gtf_select)$gene_name %in% genes_to_work_on]
raw_annotation_gtf_df_select <- as.data.frame(raw_annotation_gtf_select)
raw_annotation_gtf_df_select <- raw_annotation_gtf_df_select %>%
  dplyr::select(seqnames, start, end, width,strand,
                type, ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number) %>%
  dplyr::filter(transcript_id %in% annotation_of_genes_select$ensembl_transcript_id_version)
#### All trans GTF annotation
raw_annotation_gft_file_all_trans <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gff3"
raw_annotation_gtf_all_trans <- rtracklayer::import(con = raw_annotation_gft_file_all_trans)
raw_annotation_gtf_all_trans <- raw_annotation_gtf_all_trans[mcols(raw_annotation_gtf_all_trans)$gene_name %in% genes_to_work_on]
raw_annotation_gtf_df_all_trans <- as.data.frame(raw_annotation_gtf_all_trans)
raw_annotation_gtf_df_all_trans <- raw_annotation_gtf_df_all_trans %>%
  dplyr::select(seqnames, start, end, width,strand,
                type, ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number) %>% 
  dplyr::filter(type != "gene")

#### Remove non useful objects
remove(ensembl);remove(filts);remove(attrs)
remove(mart);remove(annotation_of_genes)
remove(raw_annotation_gtf_select);remove(raw_annotation_gtf_all_trans)

#### Add Introns
raw_annotation_gtf_df_select <- add_introns_to_gtf(gtf_to_add_info = raw_annotation_gtf_df_select)
raw_annotation_gtf_df_all_trans <- add_introns_to_gtf(gtf_to_add_info = raw_annotation_gtf_df_all_trans)

### CDS filtration
raw_annotation_gtf_df_select <- filter_CDS(
  raw_annotation = raw_annotation_gtf_df_select,
  exon_n_col = "exon_number", 
  groupping = "transcript_id")
raw_annotation_gtf_df_all_trans <- filter_CDS(
  raw_annotation = raw_annotation_gtf_df_all_trans,
  exon_n_col = "exon_number", 
  groupping = "transcript_id")

###########################################
### Flame results overlap recalculation ###
###########################################
flame_results <- flame_results %>% 
  dplyr::mutate(percent_overlap = ifelse(
    (diff_start != 0 | diff_end != 0) & percent_overlap == 1 & type %in% c("exon","CDS"),
    width/overlap_width, percent_overlap)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(id_id = paste0(raw_id,"_",id)) %>% 
  dplyr::ungroup()

## CDS filtration
flame_results <- flame_results %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(exon_number = gsub(pattern = "\\D+", replacement = "", x = raw_id)) %>% 
  dplyr::ungroup()

flame_results <- filter_CDS(
  raw_annotation = flame_results,
  exon_n_col = "exon_number", groupping = "id")

flame_results_percent_to_match <- flame_results %>% 
  dplyr::filter(info_from == "exon") %>% 
  dplyr::filter(type %in% c("CDS","exon","five_prime_UTR","three_prime_UTR")) %>% 
  dplyr::select(id_id,percent_overlap) %>% 
  dplyr::group_by(id_id) %>% 
  dplyr::mutate(check = n()) %>% 
  dplyr::ungroup()

flame_results_percent_to_match_keep <- flame_results_percent_to_match %>% 
  dplyr::group_by(id_id) %>% 
  dplyr::filter(check == 1) %>% 
  dplyr::ungroup()

flame_results_percent_to_match_remove <- flame_results_percent_to_match %>% 
  dplyr::group_by(id_id) %>% 
  dplyr::filter(check > 1) %>%
  dplyr::filter(percent_overlap != max(percent_overlap, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::pull(var = id_id)

flame_results <- flame_results %>% 
  dplyr::filter(!id_id %in% flame_results_percent_to_match_remove) %>% 
  dplyr::select(-exon_number)

remove(flame_results_percent_to_match);remove(flame_results_percent_to_match_remove)

######################################
### I Ranges GTF object generation ###
######################################

### Do the ranges
ranges_select <- ranges_generator(gff_transformed_df = raw_annotation_gtf_df_select)
# ranges_all_trans <- ranges_generator(gff_transformed_df = raw_annotation_gtf_df_all_trans)

##########################################
### Comparison with Select transcirpts ###
##########################################

### Do the comparisons
# df_comparisons_select <- comparison_maker(flame_results_df = flame_results,
#                                           ranges_to_check = ranges_select)
# 
# ### Add unknown information
# unknown_to_add <- flame_results %>% 
#   dplyr::filter(type == "Unknown")
# 
# unknown_to_add <- unknown_to_add %>% 
#   dplyr::select(overlap_width, percent_overlap, gene_symbol, transcript, sample_name, id) %>%
#   dplyr::rename(gene = gene_symbol) %>% 
#   dplyr::mutate(subject_id = "Unknown",
#                 refference_id = "Unknown",
#                 name_of_query = "Unknown") 
# 
# df_comparisons_select <- rbind(df_comparisons_select, unknown_to_add)


### Subset pefect matching Transcripts
# perfect_matching_ids <- df_comparisons_select %>% 
#   dplyr::group_by(id) %>% 
#   dplyr::filter(all(percent_overlap[!grepl(x = refference_id, pattern = "CDS")] == 1)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(id_id = paste0(refference_id,"_",id)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(percent_overlap = ifelse(id_id %in% flame_results_percent_to_match_keep$id_id,
#                                          flame_results_percent_to_match_keep$percent_overlap[match(id_id, flame_results_percent_to_match_keep$id_id)],
#                                           percent_overlap)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(id) %>% 
#   dplyr::filter(all(percent_overlap[!grepl(x = refference_id, pattern = "CDS")] == 1)) %>% 
#   dplyr::ungroup()  
  




