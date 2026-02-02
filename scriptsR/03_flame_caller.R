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
### exons data
flame_matched_exons_data <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matched_exon_data.tsv" 
### exons data
flame_matched_genes_and_transcript_data <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matched_gene_and_trans_data.tsv" 
### final_data_place
final_data_place <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/" 

######### HERE PUT THE SNAKEMAKE PARAMS SO ALL THAT IS ABOVE THIS PART SHOULD BE TAKEN BY SNAKEMAKE
######
###

#### Import the data
flame_matched_exons_data <- readr::read_tsv(file = flame_matched_exons_data)
flame_matched_genes_and_transcript_data <- readr::read_tsv(file = flame_matched_genes_and_transcript_data)

# flame_matched_exons_data <- flame_matched_exons_data %>%
#   dplyr::filter(grepl(x = transcript_id, pattern = "ENSG"))
# 
# flame_matched_genes_and_transcript_data <- flame_matched_genes_and_transcript_data %>%
#   dplyr::filter((info_from == "gene") | (grepl(x = transcript_id, pattern = "ENSG")))

#### Do the calls
final_list <- list()
for (i in 1:length(unique(flame_matched_exons_data$sample_name))) {
  #### data now
  sample_now <- unique(flame_matched_exons_data$sample_name)[i]
  flame_exons_now <- flame_matched_exons_data %>% 
    dplyr::filter(sample_name == sample_now)
  flame_genes_now <- flame_matched_genes_and_transcript_data %>% 
    dplyr::filter(sample_name == sample_now) %>% 
    dplyr::filter(info_from == "gene")
  flame_trans_now <- flame_matched_genes_and_transcript_data %>% 
    dplyr::filter(sample_name == sample_now) %>% 
    dplyr::filter(info_from == "transcript")
  
  #### work per gene
  for (g in 1:length(unique(flame_genes_now$ensembl_gene_id))) {
    gene_now <- unique(flame_genes_now$ensembl_gene_id)[g]
    flame_exons_now_gene <- flame_exons_now %>% 
      dplyr::filter(ensembl_gene_id == gene_now)
    flame_trans_now_gene <- flame_trans_now %>%
      dplyr::filter(Parent == paste0("gene:",gene_now))
    if (nrow(flame_trans_now_gene) == 0) {
      cat(paste0("MAL!","\n"))
      next
    } else if (nrow(flame_trans_now_gene) > 0) {
      for (trans in 1:length(unique(flame_trans_now_gene$raw_id))) {
        trans_now <- unique(flame_trans_now_gene$raw_id)[trans]
        flame_exons_now_trans <- flame_exons_now %>% 
          dplyr::filter(Parent == trans_now)
        
        #id level annotation
        id_of_call <- paste0(sample_now,"_",trans_now)
        trans_now <- trans_now ## aixi no m'ho deixo quan faci el df res mes!
        
        ## supporting counts annotation
        supporting_counts_of_gene <- flame_genes_now$support_count[flame_genes_now$ensembl_gene_id == gene_now]
        supporting_counts_of_transcript <- flame_trans_now$support_count[flame_trans_now$raw_id == trans_now]
        proportion_of_calls <- supporting_counts_of_transcript/supporting_counts_of_gene
        
        ## nrormal calls 
        normal_calls_exons <- flame_exons_now_trans$raw_id[flame_exons_now_trans$type != "intron" & 
                                                             flame_exons_now_trans$percent_overlap == 1]
        normal_calls_introns <- flame_exons_now_trans$raw_id[flame_exons_now_trans$type == "intron" & 
                                                               flame_exons_now_trans$percent_overlap == 0]
        normal_calls <- c(normal_calls_exons, normal_calls_introns)
        normal_calls <- paste0(x = normal_calls, collapse = "; ")
        
        ## nrormal calls 
        abnormal_calls_exons <- flame_exons_now_trans$raw_id[flame_exons_now_trans$type != "intron" &
                                                             flame_exons_now_trans$percent_overlap != 1]
        abnormal_calls_introns <- flame_exons_now_trans$raw_id[flame_exons_now_trans$type == "intron" &
                                                               flame_exons_now_trans$percent_overlap != 0]
        abnormal_calls <- c(abnormal_calls_exons, abnormal_calls_introns)
        abnormal_calls <- paste0(x = abnormal_calls, collapse = "; ")
        
        ## df of abnormality
        abnormal_calls_exons_df <- flame_exons_now_trans[(flame_exons_now_trans$type != "intron" & 
                                                               flame_exons_now_trans$percent_overlap != 1),, drop = FALSE]
        abnormal_calls_introns_df <- flame_exons_now_trans[(flame_exons_now_trans$type == "intron" & 
                                                                 flame_exons_now_trans$percent_overlap != 0),, drop = FALSE]
        
        abnormal_calls_df <- rbind(
          abnormal_calls_exons_df,
          abnormal_calls_introns_df)
        
        
        if (nrow(abnormal_calls_df) != 0) {
          abnormal_calls_df <- abnormal_calls_df %>% 
            dplyr::select(raw_id,overlap_width,percent_overlap) %>% 
            dplyr::mutate(real_lengh = round(x = (overlap_width/percent_overlap), digits = 0)) %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(final_id = paste0(raw_id,"_",overlap_width,"_",percent_overlap,"_",real_lengh)) %>% 
            dplyr::ungroup()
          abnormal_calls <- abnormal_calls_df$final_id
          abnormal_calls <- paste0(x = abnormal_calls, collapse = "; ")
        } else if (nrow(abnormal_calls_df) == 0) {
          abnormal_calls <- NA
        }
        
        ## refference trans
        reff_trans <- unique(flame_exons_now_trans$reff_transcript)
        gene_ensembl <- unique(flame_exons_now_trans$ensembl_gene_id)
        gene_symbol <- unique(flame_exons_now_trans$gene_symbol)
        
        id_of_call <- paste0(sample_now,"_",gene_symbol,"_",trans_now)
        
        called_df <- data.frame(
          "Sample" = sample_now,
          "ID_of_call" = id_of_call,
          "Flame_Transcript" = trans_now,
          "Supporting_Counts_gene" = supporting_counts_of_gene,
          "Supporting_Counts_trans" = supporting_counts_of_transcript,
          "Proportion_of_calls" = proportion_of_calls,
          "Normal_calls" = normal_calls,
          "Odd_calls" = abnormal_calls,
          "Refference_Transcript" = reff_trans,
          "Ensembl_Gene_Symbol" = gene_ensembl,
          "Gene_Symbol" = gene_symbol)
        
        id_of_call <- paste0(sample_now,"_",gene_symbol,"_",trans_now)
        
        final_list[[id_of_call]] <- called_df
        
      }
    }
  }
  }

length(final_list)
length(unique(flame_matched_exons_data$id))

final_df <- do.call(rbind, final_list)
rownames(final_df) <- NULL

rm(list = setdiff(ls(), 
                  c("final_df","final_data_place")))

openxlsx::write.xlsx(x = final_df, 
                     file = paste0(final_data_place,"calling_results.xlsx"))

#### from the final df organize the data
for (i in 1:length(unique(final_df$Sample))) {
  ## data now
  sample_now <- unique(final_df$Sample)[i]
  place_of_sample_data <- paste0(final_data_place,sample_now)
  data_now <- final_df %>% 
    dplyr::filter(Sample == sample_now)
  
  if (!dir.exists(paths = place_of_sample_data)) {
    dir.create(path = place_of_sample_data)
    }
  ## data now gene
  for (g in 1:length(unique(data_now$Gene_Symbol))) {
    gene_now <- unique(data_now$Gene_Symbol)[g]
    data_now_gene <- data_now %>% 
      dplyr::filter(Gene_Symbol == gene_now) %>% 
      dplyr::arrange(desc(Proportion_of_calls))
    openxlsx::write.xlsx(file = paste0(final_data_place,sample_now,"/",gene_now,".xlsx"), 
                         x = data_now_gene)
  }
  }



