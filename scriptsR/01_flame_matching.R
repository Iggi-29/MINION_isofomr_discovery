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
flame_base_path_results <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame" ## here they are the same but it could be different
### annotation of the genes
raw_annotation_gft_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/gff3_annotations_for_igv/gff3_filtered.gff3" ## MANE!
raw_annotation_gft_file <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/ref_files/filtered.gff3"
### genes to work on
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this.txt"

######### HERE PUT THE SNAKEMAKE PARAMS SO ALL THAT IS ABOVE THIS PART SHOULD BE TAKEN BY SNAKEMAKE
######
###

#### funnctions
### filter out not required CDSs
trashy_CDS_filtration <- function(gtf_df) {
  filtered_gft_df <- gtf_df %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::filter((exon_number == max(exon_number[type == "CDS"], na.rm = TRUE)) | 
                    (exon_number == min(exon_number[type == "CDS"], na.rm = TRUE)) | 
                    type != "CDS") %>% 
    dplyr::ungroup()
  return(filtered_gft_df)
}

### add introns
intron_addition <- function(gtf_df) {
  ## guardar el raw_gtf en una variable per a fer més fàcil en rbind
  gtf_df_to_bind <- gtf_df %>% 
    dplyr::filter(type != "transcript")
  ## primer una llista per a fer un append
  list_of_intron_additions <- list()
  ## for each transcript
  for (i in 1:length(unique(gtf_df$transcript_id))) {
    ## data now
    transcript_now <- unique(gtf_df$transcript_id)[i]
    gtf_df_now <- gtf_df %>% 
      dplyr::filter(transcript_id == transcript_now) %>% 
      dplyr::mutate(exon_number = as.numeric(exon_number)) %>% 
      dplyr::filter(type != "transcript")
    strand_now <- unique(as.character(gtf_df_now$strand))
    ## calcular en funció dels introns
    # strand (+)
    if (strand_now == "+") {
      intron_additions <- gtf_df_now %>% 
        # do the calculations
        dplyr::filter(type == "exon") %>% 
        dplyr::arrange(start) %>% 
        dplyr::mutate(
          next_exon_start = lead(start),
          intron_number = exon_number,
          intron_start = end +1,
          intron_end = next_exon_start-1) %>% 
        # filter bad introns
        dplyr::filter(!is.na(next_exon_start)) 
    } else if (strand_now == "-") {
        intron_additions <- gtf_df_now %>% 
          # do the calculations
          dplyr::filter(type == "exon") %>% 
          dplyr::arrange(desc(start)) %>% 
          dplyr::mutate(
            next_exon_end = lead(end),
            intron_number = exon_number,
            intron_start = next_exon_end +1,
            intron_end = start-1
          ) %>% 
          # filter bad introns
          dplyr::filter(!is.na(next_exon_end))
    } 
    # do the transmutation
    intron_additions <- intron_additions %>%
      dplyr::transmute(
        seqnames = seqnames,
        start = intron_start,
        end = intron_end,
        width = end - start+1,
        strand = strand,
        type = "intron",
        ID = paste0("intron:",transcript_id,":",intron_number),
        gene_id = gene_id,
        gene_name = gene_name,
        transcript_id = transcript_id,
        exon_id = exon_id,
        exon_number = intron_number)
    list_of_intron_additions[[i]] <- intron_additions
  }
  ## fem la llista
  list_of_intron_additions <- do.call(rbind, list_of_intron_additions)
  ## fem el rbind amb l'original i sort per start
  gtf_w_introns <- rbind(list_of_intron_additions, gtf_df_to_bind)
  gtf_w_introns <- gtf_w_introns %>%
    dplyr::arrange(gene_name,transcript_id,start) %>% 
    dplyr::mutate(exon_number = as.numeric(exon_number))
  return(gtf_w_introns)
  }

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

##### BiomaRt data importation 
#### Connect to the service
mart <- biomaRt::useEnsembl(biomart = "ensembl", version = 115)
datasets <- biomaRt::listDatasets(mart = mart)
ensembl <- biomaRt::useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
#### Get the filters and attributes
filts <- biomaRt::listFilters(mart = ensembl)
attrs <- biomaRt::listAttributes(mart = ensembl)
#### Annotations of genes
### Annotation of genes
annotation_of_genes <- biomaRt::getBM(
  mart = ensembl,
  attributes = c("ensembl_gene_id","ensembl_gene_id_version",
                 "hgnc_symbol","ensembl_transcript_id_version"),
  filters = c("hgnc_symbol"),
  values = genes_to_work_on)
### Annotation of genes Select 
annotation_of_genes_select <- biomaRt::getBM(
  mart = ensembl,
  attributes = c("hgnc_symbol","ensembl_transcript_id_version","transcript_mane_select"),
  filters = c("hgnc_symbol"), 
  values = genes_to_work_on)
annotation_of_genes_select <- annotation_of_genes_select %>% 
  dplyr::filter(transcript_mane_select != "")
annotation_of_genes <- annotation_of_genes %>% 
  dplyr::mutate(is_MANE_select = ifelse(
    ensembl_transcript_id_version %in% annotation_of_genes_select$ensembl_transcript_id_version,
    "yes","no"))

remove(annotation_of_genes_select);remove(attrs);remove(datasets);remove(filts)
remove(mart);remove(ensembl)

##### GTF annotation work
#### GTF all the information
raw_annotation_gtf <- rtracklayer::import(con = raw_annotation_gft_file)
raw_annotation_gtf <- raw_annotation_gtf[mcols(raw_annotation_gtf)$transcript_id %in% annotation_of_genes$ensembl_transcript_id_version]
raw_annotation_gtf <- as.data.frame(raw_annotation_gtf) %>% 
  dplyr::select(seqnames,start,end,width,strand,
                type,ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number)

all(unique(raw_annotation_gtf$transcript_id) %in% annotation_of_genes$ensembl_transcript_id_version)

# #### GTF strand +
# gtf_MSH6 <- raw_annotation_gtf[raw_annotation_gtf$gene_name == "MSH6",, drop = FALSE]
# gtf_MSH6 <- gtf_MSH6[gtf_MSH6$type != "transcript",, drop = FALSE]
# gtf_MSH6$exon_number <- as.numeric(gtf_MSH6$exon_number)
# 
# gtf_MSH6 <- trashy_CDS_filtration(gtf_df = gtf_MSH6)
# gtf_MSH6 <- intron_addition(gtf_df = gtf_MSH6)
# 
# #### GTF strand -
# gtf_PALB2 <- raw_annotation_gtf[raw_annotation_gtf$gene_name == "PALB2",, drop = FALSE]
# gtf_PALB2 <- gtf_PALB2[gtf_PALB2$type != "transcript",, drop = FALSE]
# gtf_PALB2$exon_number <- as.numeric(gtf_PALB2$exon_number)
# 
# gtf_PALB2 <- trashy_CDS_filtration(gtf_df = gtf_PALB2)
# gtf_PALB2 <- intron_addition(gtf_df = gtf_PALB2)

#### with the whole GTF
raw_annotation_gtf <- trashy_CDS_filtration(gtf_df = raw_annotation_gtf)
raw_annotation_gtf <- intron_addition(gtf_df = raw_annotation_gtf)
raw_annotation_gtf <- raw_annotation_gtf %>% 
  dplyr::mutate(is_MANE_select = 
                  ifelse(transcript_id %in% annotation_of_genes$ensembl_transcript_id_version[annotation_of_genes$is_MANE_select == "yes"],
                         "yes","no"))

#### Only gene and transcript annotation
# raw_annotation_gft_genes <- rtracklayer::import(con = raw_annotation_gft_file)
# raw_annotation_gft_genes <- raw_annotation_gft_genes[mcols(raw_annotation_gft_genes)$gene_name %in% genes_to_work_on]
# raw_annotation_gft_genes <- as.data.frame(raw_annotation_gft_genes) %>%
#   dplyr::select(seqnames,start,end,width,strand,
#                 type, ID,
#                 gene_id,gene_name,
#                 transcript_id,
#                 exon_id,exon_number) %>%
#   dplyr::select(ID, gene_id, gene_name, start, end, type)
# 
# raw_annotation_gtf_transcripts <- raw_annotation_gft_genes %>% 
#   dplyr::filter(type == "transcript")
# raw_annotation_gft_genes <- raw_annotation_gft_genes %>% 
#   dplyr::filter(type == "gene")

##### FLAME output work
matching_all <- list()
for (i in 1:length(flame_results)) {
  #### Data now
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
  cat(paste0("Working on sample ",sample_name," which is sample ",i," out of ",length(flame_results),"\n"))
  
  # flame_results_now <- flame_results_now %>% 
  #   dplyr::filter(gene_id == "ENSG00000132781.22" | 
  #                   transcript_id == "ENST00000456914.7" | 
  #                   grepl(x = Parent, pattern = "ENST00000456914.7"))
  
  ### Get the genes with at least a call in our dataset
  flame_genes <- flame_results_now$ID[!is.na(flame_results_now$gene_id)]
  flame_genes <- sapply(X = flame_genes, FUN = function(x) {gsub(pattern = "gene\\:",replacement = "", x = x)},
                        USE.NAMES = FALSE)
  # flame_genes <- sapply(X = flame_genes, FUN = function(x) {gsub(pattern = "\\..*$", replacement = "", x = x)},
  #                       USE.NAMES = FALSE)
  ### Work for each of the genes
  for (g in 1:length(flame_genes)) {
    flame_gene_now <- flame_genes[g]
    gene_symbol_now <- unique(annotation_of_genes$hgnc_symbol[match(
      gsub(pattern = "\\..*$", replacement = "", x = flame_gene_now), 
      annotation_of_genes$ensembl_gene_id)])
    cat(paste0("Working on gene ",gene_symbol_now," which is gene ",g," out of ",length(flame_genes),"\n"))
    
    ### Pick the data for that gene <- FLAME level
    flame_results_now_gene <- flame_results_now %>% 
      dplyr::filter(type == "gene") %>% 
      dplyr::filter(gene_id == flame_gene_now)
    flame_results_now_trans <- flame_results_now %>% 
      dplyr::filter(type == "transcript") %>% 
      dplyr::filter(grepl(x = Parent,pattern = flame_gene_now)) %>% 
      dplyr::pull(ID)
    
    ### Pick the data for that gene <- Raw GTF annotation data
    raw_annotation_gtf_now <- raw_annotation_gtf %>% 
      dplyr::filter(gene_name == gene_symbol_now) %>% 
      dplyr::filter(type != "transcript")
    
    ### Work the transcripts
    for (trans in 1:length(flame_results_now_trans)) {
      ### Data now
      transcript_now <- flame_results_now_trans[trans]
      ## only transcript data
      flame_transcript_now_data <- flame_results_now %>% 
        dplyr::filter(type == "transcript") %>% 
        dplyr::filter(ID == transcript_now)
      ## only exon data of that transcript
      flame_exon_now_data <- flame_results_now %>% 
        dplyr::filter(Parent == transcript_now) 
      ## flame exon ranges
      exon_ranges <- IRanges::IRanges(start = flame_exon_now_data$start,
                                      end = flame_exon_now_data$end,
                                      names = flame_exon_now_data$exon_id)
      
      ## Match each of the possible transcripts of the ref with the 
      for (reff_trans in 1:length(unique(raw_annotation_gtf_now$transcript_id))) {
        ### Data now
        reff_trans_now <- unique(raw_annotation_gtf_now$transcript_id)[reff_trans]
        raw_annotation_gtf_now_trans <- raw_annotation_gtf_now %>% 
          dplyr::filter(transcript_id == reff_trans_now)  
        raw_ranges <- IRanges::IRanges(start = raw_annotation_gtf_now_trans$start,
                                       end = raw_annotation_gtf_now_trans$end,
                                       names = raw_annotation_gtf_now_trans$ID)
        
        ### Do the IRanges opperation
        hits <- IRanges::findOverlaps(subject = exon_ranges,
                                      query = raw_ranges,
                                      minoverlap = 0)
        
        if (length(hits) > 0) {
          overlaps <- pintersect(x = exon_ranges[subjectHits(hits)],
                                 y = raw_ranges[queryHits(hits)])
          percent_overlap <- width(overlaps)/width(raw_ranges[queryHits(hits)])
          
          ## matching df
          matching_df <- data.frame(
            raw_id = names(raw_ranges)[queryHits(hits)],
            exon_id = names(exon_ranges)[subjectHits(hits)],
            overlap_width = width(overlaps),
            percent_overlap = round(x = percent_overlap, digits = 4))
          
          ## add non existing matches at the refference
          non_matched_refference <- setdiff(raw_annotation_gtf_now_trans$ID,
                                 matching_df$raw_id)
          if (length(non_matched_refference) > 0) {
            non_matched_df_refference <- data.frame(
              raw_id = non_matched_refference,
              exon_id = "NO PRESENT",
              overlap_width = 0,
              percent_overlap = 0)
            matching_df <- rbind(matching_df,non_matched_df_refference)
          }
          
          ## add non existing matches at the flame_result
          non_matched_flame <- setdiff(flame_exon_now_data$exon_id,
                                            matching_df$exon_id[matching_df$exon_id != "NO PRESENT"])
          if (length(non_matched_flame) > 0) {
            non_matched_df_flame <- data.frame(
              raw_id = "NO PRESENT",
              exon_id = non_matched_flame,
              overlap_width = 0,
              percent_overlap = 0)
            matching_df <- rbind(matching_df,non_matched_df_flame)
          }
          
          ## Make the matching_df a trackable data.frame
          matching_df$flame_transcript <- transcript_now
          matching_df$reff_trans_now <- reff_trans_now
          matching_df$sample_name <- sample_name
          id_of_comparison <- paste0(transcript_now,"_",reff_trans_now,"_",sample_name)
          matching_df$id_of_comparison <- id_of_comparison
          matching_all[[id_of_comparison]] <- matching_df }
        else if (length(hits)) {
          next}
        }
      }
    }
  }

#### End all the matching data
matching_all_df <- do.call(rbind, matching_all)
rownames(matching_all_df) <- NULL

write.table(x = matching_all_df, 
            file = "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/raw_matched_df.txt", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

#### Select the best match per flameID
matching_all_df <- matching_all_df %>% 
  ## reasign type of coincidence
  dplyr::rowwise() %>% 
  dplyr::mutate(type = gsub(pattern = "\\:.*$", replacement = "", x = raw_id)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(type = ifelse(raw_id == "NO PRESENT",
                              "extent",type)) %>% 
  dplyr::relocate(type, .after = exon_id) %>% 
  ## do the lenght of extent
  dplyr::mutate(length_of_ext = ifelse(type == "extent",
                                      abs(as.numeric(strsplit(x = gsub(x = exon_id, pattern = "exon\\:", replacement = ""), split = "_")[[1]][1]) - 
                                            as.numeric(strsplit(x = gsub(x = exon_id, pattern = "exon\\:", replacement = ""), split = "_")[[1]][2])),
                                      0)) %>% 
  ## mean overlap 
  dplyr::group_by(sample_name,flame_transcript,reff_trans_now) %>%
  dplyr::mutate(mean_overlap_non_int = mean(percent_overlap[!type %in% c("intron","extent")]),
                mean_overlap_int = mean(percent_overlap[type %in% "intron"]),
                mean_length_of_ext = mean(length_of_ext)) %>% 
  dplyr::ungroup()

### filter best match per comparison
matching_all_df_filterer <- matching_all_df %>%
  dplyr::select(
    sample_name, flame_transcript, reff_trans_now,
    mean_overlap_non_int,mean_overlap_int,mean_length_of_ext,
    id_of_comparison) %>%
  dplyr::distinct() %>% 
  ## do the filtration
  dplyr::group_by(sample_name,flame_transcript) %>% 
  dplyr::arrange(
    desc(mean_overlap_non_int),
    mean_overlap_int,
    mean_length_of_ext,
    .by_group = TRUE) %>%
  dplyr::slice(1) %>% 
  dplyr::ungroup()

write.table(x = matching_all_df_filterer, 
            file = "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matching_scores.txt", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### matching df final
matching_final_df <- matching_all_df %>% 
  # filter the best comparison only and remove redundant columns
  dplyr::filter(id_of_comparison %in% matching_all_df_filterer$id_of_comparison) %>% 
  dplyr::select(sample_name,exon_id,raw_id,type,overlap_width,percent_overlap,flame_transcript,reff_trans_now) %>% 
  dplyr::rename(reff_transcript = reff_trans_now)

write.table(x = matching_final_df, 
            file = "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/matching_to_annotate.txt", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

openxlsx::write.xlsx(x = raw_annotation_gtf,
                     file = "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame/final/annotation_of_transcripts.xlsx")






