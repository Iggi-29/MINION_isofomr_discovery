#######################################
### FLAMES long reads data analysis ###
#######################################

#### libraries
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

place_to_save_the_data <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flames"
if(!dir.exists(paths = paste0(place_to_save_the_data,"/final"))) {
  dir.create(path = paste0(place_to_save_the_data,"/final"))
  cat(paste0("Created the final data dir","\n"))
} else {
  cat(paste0("Didn't need to create the final data dir","\n"))
}

#### functions
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

#### annotation data
raw_annotation_gft_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/gff3_annotations_for_igv/gff3_filtered.gff3"
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this_ENIGMA.txt"

#### flame results
flame_base_path_results <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flames"
flame_results <- list.files(path = flame_base_path_results, 
                            recursive = TRUE, full.names = TRUE,
                            pattern = "isoform_annotated\\.filtered\\.gff3$")
flame_results <- flame_results[grep(x = flame_results, pattern = "LRS4_LRS4_barcode01|LRS4_LRS4_barcode02")]

#### Biomart data importation
genes_to_work_on <- readLines(con = genes_to_work_on)
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
                                                     "hgnc_symbol","ensembl_transcript_id_version"),
                                      filters = c("hgnc_symbol"), 
                                      values = genes_to_work_on)
annotation_of_genes_select <- biomaRt::getBM(mart = ensembl,
                                             attributes = c("hgnc_symbol","ensembl_transcript_id_version","transcript_mane_select"), 
                                             filters = c("hgnc_symbol"), values = genes_to_work_on)
annotation_of_genes_select <- annotation_of_genes_select %>% 
  dplyr::filter(transcript_mane_select != "")

annotation_of_genes <- annotation_of_genes %>% 
  dplyr::mutate(is_MANE_select = ifelse(ensembl_transcript_id_version %in% annotation_of_genes_select$ensembl_transcript_id_version,
                                        "yes","no"))

#### gtf annotation work
### genes
raw_annotation_gtf_genes <- rtracklayer::import(con = raw_annotation_gft_file)
raw_annotation_gtf_genes <- raw_annotation_gtf_genes[mcols(raw_annotation_gtf_genes)$gene_name %in% genes_to_work_on]
raw_annotation_gtf_genes <- as.data.frame(raw_annotation_gtf_genes) %>% 
  dplyr::select(seqnames, start, end, width,strand, 
                type, ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number) %>% 
  dplyr::select(ID,gene_id,gene_name,start,end,type)

### transcripts and exons
raw_annotation_gtf <- rtracklayer::import(con = raw_annotation_gft_file)
raw_annotation_gtf <- raw_annotation_gtf[mcols(raw_annotation_gtf)$transcript_id %in% annotation_of_genes$ensembl_transcript_id_version]
raw_annotation_gtf_df <- as.data.frame(raw_annotation_gtf)
raw_annotation_gtf_df <- raw_annotation_gtf_df %>% 
  dplyr::select(seqnames, start, end, width,strand, 
                type, ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number)

all(unique(raw_annotation_gtf_df$transcript_id) %in% annotation_of_genes$ensembl_transcript_id_version)

### transcripts only
raw_annotation_gtf_transcripts <- raw_annotation_gtf_genes %>% 
  dplyr::filter(type == "transcript") 

### Add intron information
raw_annotation_gtf_df <- add_introns_to_gtf(gtf_to_add_info = raw_annotation_gtf_df)

all_data_list <- list()
for (i in 1:length(flame_results)) {
  #### Data_now
  flame_results_now_filename <- flame_results[i]
  flame_results_now <- rtracklayer::import(con = flame_results_now_filename)
  flame_results_now_df <- as.data.frame(flame_results_now)
  flame_results_now_df$Parent <- sapply(
    flame_results_now_df$Parent,
    function(x) if (length(x) == 0) NA_character_ else x[1]
  )  

  sample_name <- gsub(pattern = ".*/(.*)/[^/]+$", 
                      replacement = "\\1", 
                      x = flame_results_now_filename)
  cat(paste0("Working on sample ",sample_name," which is sample ",i," out of ",length(flame_results),"\n"))
  
  #### Genes of data now
  flame_genes <- flame_results_now_df$ID[!is.na(flame_results_now_df$gene_id)]
  flame_genes <- sapply(X = flame_genes, FUN = function(x) {gsub(pattern = "gene\\:",replacement = "", x = x)},
                        USE.NAMES = FALSE)
  flame_genes <- sapply(X = flame_genes, FUN = function(x) {gsub(pattern = "\\..*$", replacement = "", x = x)},
                        USE.NAMES = FALSE)
  #### Work for each of the genes
  for (g in 1:length(flame_genes)) {
    ##### Data now
    flame_genes_now <- flame_genes[g]
    gene_symbol <- unique(annotation_of_genes$hgnc_symbol[match(flame_genes_now, annotation_of_genes$ensembl_gene_id)])
    #### Reconstruct the data
    gene_data <- flame_results_now_df %>%
      dplyr::filter(type == "gene") %>%
      dplyr::filter(grepl(x = ID, pattern = flame_genes_now))
    transcript_data <- flame_results_now_df %>%
      dplyr::filter(type == "transcript") %>%
      dplyr::filter(grepl(x = Parent, pattern = flame_genes_now))
    
    #### Work for each of the transcripts
    for (t in 1:length(transcript_data$ID)) {
      ##### Data now
      transcript_now <- transcript_data$ID[t]
      transcript_data_now <- transcript_data %>%
        dplyr::filter(ID == transcript_now)
      #### Reconstruct the data
      exon_data <- flame_results_now_df %>%
        dplyr::filter(type == "exon") %>%
        dplyr::filter(Parent == transcript_now) %>% 
        dplyr::mutate(rank = as.numeric(rank)) %>%
        dplyr::select(-c(type))

      ##### Which exons are the exons ?
      #### work the raw_gtf_data
      raw_annotation_gtf_df_now <- raw_annotation_gtf_df %>%
        dplyr::filter(gene_name == gene_symbol) %>%
        dplyr::filter(type != "transcript") %>%
        dplyr::filter(transcript_id == annotation_of_genes_select$ensembl_transcript_id_version[annotation_of_genes_select$hgnc_symbol == gene_symbol]) %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(id_of_app = paste0(type,"_",exon_number)) %>%
        dplyr::ungroup() %>%
        dplyr::select(id_of_app,type,start,end)

      #### Generate IRanges data
      exon_ranges <- IRanges(start = exon_data$start,
                             end = exon_data$end,
                             names = exon_data$exon_id)

      raw_ranges <- IRanges(start = raw_annotation_gtf_df_now$start,
                            end = raw_annotation_gtf_df_now$end,
                            names = raw_annotation_gtf_df_now$id_of_app)

      #### Calculate the ranges
      hits <- findOverlaps(subject = exon_ranges,
                           query = raw_ranges,
                           minoverlap = 0)
      overlaps <- pintersect(x = exon_ranges[subjectHits(hits)],
                             y = raw_ranges[queryHits(hits)])
      percent_overlap <- width(overlaps)/width(raw_ranges[queryHits(hits)])

      #### Matching_df
      matching_df <- data.frame(
        raw_id = names(raw_ranges)[queryHits(hits)],
        exon_id = names(exon_ranges)[subjectHits(hits)],
        overlap_width = width(overlaps),
        percent_overlap = round(x = percent_overlap, digits = 4))

      ### Get original coordinates for Matching df
      raw_annotation_gtf_df_now <- raw_annotation_gtf_df_now %>%
        dplyr::rename(raw_id = id_of_app)
      matching_df <- matching_df %>%
        dplyr::left_join(y = raw_annotation_gtf_df_now, by = "raw_id") %>%
        dplyr::rename(original_start = start,
                      original_end = end)
      
      exon_data <- dplyr::left_join(x = exon_data, y = matching_df, by = "exon_id")
      exon_data <- exon_data %>%
        dplyr::relocate(type, .after = source) %>%
        dplyr::mutate(type = as.character(type)) %>% 
        dplyr::mutate(
          raw_id = ifelse(is.na(raw_id),"Unknown",raw_id),
          type = ifelse(is.na(type),"Unknown",type),
          overlap_width = ifelse(is.na(overlap_width),0,overlap_width),
          percent_overlap = ifelse(is.na(percent_overlap),0,percent_overlap)) %>% 
        dplyr::ungroup() %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(raw_id = ifelse(raw_id == "Unknown", paste0(x = c(start,"_",raw_id,"_",end), collapse = ""),raw_id)) %>%
        dplyr::ungroup()

      ### Finish the all_data df
      ## Gene data
      original_gene <- raw_annotation_gtf_genes %>%
        dplyr::filter(type == "gene") %>% 
        dplyr::filter(gene_id == unique(gene_data$gene_id)) %>%
        dplyr::rename(original_start = start,
                      original_end = end)
      
      gene_data <- gene_data %>%
        dplyr::mutate(overlap_width = NA,
                      percent_overlap = NA,
                      raw_id = "",
                      original_start = unique(original_gene$original_start),
                      original_end = unique(original_gene$original_end))
      
      ## Transcript data
      original_transcript_id <- unique(transcript_data_now$ID)
      original_transcript_id <- gsub(pattern = "transcript\\:", replacement = "", x = original_transcript_id)
      original_transcript_id <- gsub(pattern = "_.*", replacement = "", x = original_transcript_id)
      
      if (grepl(x = original_transcript_id, pattern = "ENSG")) {
        original_transcript <- raw_annotation_gtf_transcripts %>%
          dplyr::filter(gene_id == original_transcript_id) %>%
          dplyr::filter(ID == annotation_of_genes_select$ensembl_transcript_id_version[annotation_of_genes_select$hgnc_symbol == gene_symbol]) %>% 
          dplyr::rename(original_start = start,
                        original_end = end)  

      } else if (grepl(x = original_transcript_id, pattern = "ENST")) {
        original_transcript <- raw_annotation_gtf_transcripts %>%
          dplyr::filter(ID == original_transcript_id) %>%
          dplyr::rename(original_start = start,
                        original_end = end)}

      transcript_data_now <- transcript_data_now %>%
        dplyr::mutate(overlap_width = NA,
                      percent_overlap = NA,
                      raw_id = "",
                      original_start = unique(original_transcript$original_start),
                      original_end = unique(original_transcript$original_end))


      gene_data <- gene_data %>% mutate(rank = as.character(rank))
      transcript_data_now <- transcript_data_now %>% mutate(rank = as.character(rank))
      exon_data <- exon_data %>% mutate(rank = as.character(rank))

      gene_data$info_from <- "gene"
      transcript_data_now$info_from <- "trans"
      exon_data$info_from <- "exon"

      all_data <- bind_rows(gene_data, transcript_data_now, exon_data)
      all_data <- all_data %>%
        dplyr::select(-c(score, phase)) %>%
        dplyr::mutate(ID = ifelse(is.na(ID), exon_id,ID)) %>%
        dplyr::select(-c(exon_id,transcript_id,exon_id))

      ### All data concatenation
      unique_id <- paste0(sample_name,"_",gene_symbol,"_",transcript_now)
      all_data <- all_data %>%
        dplyr::mutate(id = unique_id,
                      gene_symbol = gene_symbol,
                      transcript = transcript_now,
                      sample_name = sample_name) %>%
        dplyr::relocate(gene_symbol, transcript, .after = gene_id) %>%
        dplyr::relocate(sample_name, .before = 1)

      exon_range <- all_data$raw_id
      exon_range <- exon_range[grep(x = exon_range, pattern = "exon")]
      exon_range <- sapply(X = exon_range, FUN = function(x) {
        exon_num <- gsub(pattern = "exon_", replacement = "", x = x)
        exon_num <- as.numeric(exon_num)
        return(exon_num)
      }, USE.NAMES = FALSE)
      exon_range <- range(exon_range)
      exon_range <- paste0(x = c(as.character(exon_range)), collapse = "-")

      all_data <- all_data %>%
        dplyr::mutate(exon_range = exon_range) %>%
        dplyr::relocate(exon_range, .after = support_count) %>%
        dplyr::relocate(original_start, original_end, .after = end)
      
      all_data_list[[unique_id]] <- all_data

      }
    }
  }

#### Have all the dfs finished and save the data
all_data_def <- do.call(rbind, all_data_list)
rownames(all_data_def) <- NULL

all_data_def <- all_data_def %>%
  dplyr::rowwise() %>%
  dplyr::mutate(original_start = ifelse(type == "gene", raw_annotation_gtf_genes$start[match(gene_id, raw_annotation_gtf_genes$gene_id)], original_start)) %>%
  dplyr::mutate(original_end = ifelse(type == "gene", raw_annotation_gtf_genes$end[match(gene_id, raw_annotation_gtf_genes$gene_id)], original_end)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(diff_start = start - original_start,
                diff_end = end - original_end) %>%
  dplyr::relocate(diff_start, diff_end, .after = original_end)

#### Overlapping correction
# all_data_def <- all_data_def %>%
#   dplyr::mutate(percent_overlap = ifelse(
#     (diff_start != 0 | diff_end != 0) & percent_overlap == 1 & type %in% c("exon","CDS",
#                                                                            "five_prime_UTR","three_prime_UTR"),
#     width/overlap_width, percent_overlap))  

#### Garbage CDS removal
all_data_def <- all_data_def %>%
  dplyr::mutate(rank = as.numeric(rank)) %>%
  dplyr::group_by(id) %>%
  dplyr::filter(
    type != "CDS" |
      rank == dplyr::if_else(
        any(type == "CDS"),
        max(rank[type == "CDS"], na.rm = TRUE),
        NA_real_
      ) |
      rank == dplyr::if_else(
       any(type == "CDS"),
       min(rank[type == "CDS"], na.rm = TRUE),
       NA_real_)
  ) %>% 
  dplyr::ungroup()
  
openxlsx::write.xlsx(
  x = all_data_def,
  file = paste0(place_to_save_the_data,"/final/annotated_flames_data2.xlsx"))
readr::write_tsv(
  x = all_data_def,
  file = paste0(place_to_save_the_data,"/final/annotated_flames_data2.tsv"))

