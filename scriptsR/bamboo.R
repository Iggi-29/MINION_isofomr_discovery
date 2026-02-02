######################################
### Bambu long reads data analysis ###
######################################

#### libraries
library(bambu)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

subdirs_to_generate <- c("gtf","gtf_canonical")

#### Run depending information
### Place to save the data
place_to_save_the_data <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/bambu"
# place_to_save_the_data <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run1/results/bambu/"
### Genes to work on
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this_ENIGMA.txt"
# genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this.txt"
### Alignment data
base_path_alignments <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/alignments/flair"
# base_path_alignments <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run1/results/alignments/flair/"
### List of bathces
list_of_batches <- c("LRS4_LRS4","LRS5_LRS5","LRS6_LRS6","LRS7_LRS7")
# list_of_batches <- c("Run_17112025")

#### Generate the folders just in case
for (i in 1:length(subdirs_to_generate)) {
  dir_now <- subdirs_to_generate[i]
  if (!dir.exists(paths = paste0(place_to_save_the_data,"/",dir_now))) {
    dir.create(path = paste0(place_to_save_the_data,"/",dir_now))
    cat(paste0(dir_now," generated","\n"))
  } else {
    cat(paste0(dir_now," did already exist","\n"))
  }
}


#### Pick the raw fasta file and gtf 
fa_file <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/GRCh38.primary_assembly.genome.fa"
gtf_file_name <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf"
#### Pick the alignment files
list_of_files <- list.files(path = base_path_alignments, pattern = "\\.flair\\.aligned\\.bam$", full.names = TRUE)

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


#### pick the gtf and filter it 
gtf_file_raw <- rtracklayer::import(con = gtf_file_name)
gtf_file_raw <- gtf_file_raw[gtf_file_raw$type == "CDS"]
gtf_file_raw <- gtf_file_raw[
  mcols(gtf_file_raw)$gene_name %in%
    annotation_of_genes$hgnc_symbol]

annotation_of_genes <- annotation_of_genes %>% 
  dplyr::filter(ensembl_transcript_id_version %in% gtf_file_raw$transcript_id)

for (batch_now in 1:length(list_of_batches)) {
  ##### Data now
  batch_to_work <- list_of_batches[batch_now]
  bam_files_now <- grep(pattern = batch_to_work, x = list_of_files, value = TRUE)
  cat(paste0("Working on batch: ",batch_to_work," which is ",batch_now," out of ",length(list_of_batches),"\n"))
  
  #### Generate the Bambu data
  # gtf file
  tictoc::tic()
  gtf_file <- bambu::prepareAnnotations(x = gtf_file_name)
  tictoc::toc()
  # bambu data w the bams
  tictoc::tic()
  se <- bambu::bambu(reads = bam_files_now, annotations = gtf_file, genome = fa_file, NDR = 1)
  tictoc::toc()
  
  #### Check the Bambu data
  # check the quantification levels
  assayNames(se)
  # check the metadata of our exeriment
  head(colData(se))
  # The annotation of each found transcript 
  head(mcols(rowRanges(se)))
  # Check the number of transcripts
  table(mcols(rowRanges(se))$GENEID) |> summary()
  
  #### Filter junk isoforms
  ## Extract assays
  counts <- assay(se, "counts")
  flc <- assay(se, "fullLengthCounts")
  
  ## Filtering rule
  keep_tx <- (rowSums(flc > 0) >= 1) &      # full-length support in ≥1 samples
    (rowSums(flc) >= 3) &          # ≥3 full-length reads in total
    (rowSums(counts) >= 25)        # sufficient overall expression
  ## Apply filter
  se_robust <- se
  # se_robust <- se[keep_tx, ]
  
  #### Compare isoforms
  rr <- rowRanges(se_robust)
  #### Genes with canonical isoforms
  annotated_isoforms <- as.data.frame(rr@elementMetadata)
  
  #### Filter the isoforms data
  annotated_isoforms <- annotated_isoforms %>% 
    ## filter for non-novel genes
    dplyr::filter(novelGene == FALSE) %>% 
    ## add check if isoforms are in the genes of interest and filter
    dplyr::mutate(normal_gene_id = gsub(pattern = "\\..*$", replacement = "", x = GENEID)) %>% 
    dplyr::mutate(Gene_in_interest = normal_gene_id %in% annotation_of_genes$ensembl_gene_id) %>% 
    dplyr::relocate(Gene_in_interest, normal_gene_id, .before = GENEID) %>% 
    dplyr::rename(GENEID_ensembl = normal_gene_id) %>% 
    dplyr::filter(Gene_in_interest == TRUE) %>% 
    ## add the annotations of the genes 
    dplyr::mutate(Gene_symbol = annotation_of_genes$hgnc_symbol[match(GENEID_ensembl,
                                                                      annotation_of_genes$ensembl_gene_id)]) %>% 
    dplyr::relocate(Gene_symbol)
  
  final_data_list <- list()
  gtf_by_tx <- split(gtf_file_raw, mcols(gtf_file_raw)$transcript_id)
  
  #### For each gene, do the check
  for (i in 1:length(unique(annotated_isoforms$Gene_symbol))) {
    ### data_now
    gene_now <- unique(annotated_isoforms$Gene_symbol)[i]
    data_now <- annotated_isoforms %>% 
      dplyr::filter(Gene_symbol == gene_now)
    if (all(data_now$TXNAME %in% gtf_file_raw$transcript_id)) {
      cat(paste0("All isoforms of ",gene_now," are canonical No comparison will be performed","\n"))
    } else if (!all(data_now$TXNAME %in% gtf_file_raw$transcript_id)) {
      cat(paste0("Not all isoforms of ",gene_now," are canonical, comparsion will be performed","\n"))
      # pick reference GRAnges to compare
      transcripts_to_compare_with <- annotation_of_genes %>%
        dplyr::filter(hgnc_symbol == gene_now) %>%
        dplyr::pull(ensembl_transcript_id_version)
      
      ranges_of_subject <- rr[mcols(rr)$TXNAME %in% transcripts_to_compare_with]
      ranges_of_query <- rr[mcols(rr)$TXNAME %in% data_now$TXNAME[!data_now$TXNAME %in% transcripts_to_compare_with]]
      
      for (j in 1:length(names(ranges_of_query))) {
        ## data now
        query_now <- names(ranges_of_query)[j]
        ranges_of_query_now <- ranges_of_query[mcols(ranges_of_query)$TXNAME %in% query_now]
        for (k in 1:length(names(ranges_of_subject))) {
          ## data now
          subject_now <- names(ranges_of_subject)[k]
          ranges_of_subject_now <- ranges_of_subject[mcols(ranges_of_subject)$TXNAME %in% subject_now]
          ## final data
          name_final <- paste0(query_now,"_",subject_now)
          comparison_now <- bambu::compareTranscripts(query = ranges_of_query_now,
                                                      subject = ranges_of_subject_now)
          final_data_list[[name_final]] <- comparison_now
        }
      }
    }
  }
  
  final_data <- do.call(rbind, final_data_list)
  final_data <- merge(x = final_data, y = annotation_of_genes, 
                      by.x = "subjectId", by.y = "ensembl_transcript_id_version", 
                      all.x = TRUE)
  final_data <- final_data %>% 
    dplyr::relocate(ensembl_gene_id_version, hgnc_symbol) %>% 
    dplyr::select(-c(ensembl_gene_id)) %>% 
    dplyr::relocate(is_MANE_select, .after = subjectId)
  
  ##### get the expression
  final_data_ids <- unique(x = c(final_data$queryId, final_data$subjectId))
  filtered_counts <- counts[grep(pattern = paste0(x = final_data_ids, collapse = "|"), x = rownames(counts)),,
                            drop = FALSE]
  
  ##### get the data to visualize
  se_export <- se_robust[mcols(rowRanges(se_robust))$TXNAME %in% final_data_ids,]
  rr_export <- rowRanges(se_export)
  
  se_export_canonical <- se_export[
    mcols(rowRanges(se_export))$TXNAME %in% c(
      final_data_ids[final_data_ids %in% annotation_of_genes_select$ensembl_transcript_id_version],
      final_data_ids[!final_data_ids %in% annotation_of_genes$ensembl_transcript_id_version])]
  rr_export_canonical <- rowRanges(se_export_canonical)
  
  #### final counts
  final_counts <- as.data.frame(filtered_counts)
  final_counts$queryId <- rownames(final_counts)
  final_counts <- final_counts %>% 
    dplyr::relocate(queryId, .before = 1)
  rownames(final_counts) <- NULL
  final_counts <- final_counts %>% 
    dplyr::filter(queryId %in% final_data$queryId)
  colnames(final_counts)[2:ncol(final_counts)] <- sapply(
    X = colnames(final_counts)[2:ncol(final_counts)],
    FUN = function(x) {
      gsub(
        pattern = "\\.flair\\.aligned",
        replacement = "", x = x)})
  
  final_data2 <- dplyr::left_join(x = final_data, y = final_counts, by = "queryId")
  
  ##### save the data##### sNULLave the data
  bambu::writeToGTF(annotation = rr_export,
                    file = paste0(place_to_save_the_data,"/gtf/",batch_to_work,".gtf"))
  
  bambu::writeToGTF(annotation = rr_export_canonical,
                    file = paste0(place_to_save_the_data,"/gtf_canonical/",batch_to_work,".gtf"))
  
  
  openxlsx::write.xlsx(x = final_data2, 
                       file = paste0(place_to_save_the_data,"/",batch_to_work,"_results",".xlsx"))
  }




