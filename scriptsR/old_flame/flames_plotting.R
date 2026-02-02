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

#### plots dir
flame_plots_dir <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flames/final/plots"
if (!dir.exists(paths = flame_plots_dir)) {
  dir.create(flame_plots_dir)
  cat(paste0("The plots dir has been created. \n"))
} else {
  cat(paste0("The plots dir already exists, no further action has been taken. \n"))
}

#### flame results
flame_results_file_name <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flames/final/annotated_flames_data2.xlsx"
flame_results <- openxlsx::read.xlsx(xlsxFile = flame_results_file_name)
flame_results <- flame_results %>% 
  dplyr::filter(gene_symbol == "MUTYH") %>% 
  dplyr::filter(sample_name == "LRS4_LRS4_barcode01")

#### genes to work on
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this_ENIGMA.txt"

#### Biomart data importation
genes_to_work_on <- readLines(con = genes_to_work_on)
genes_to_work_on <- "MUTYH"
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

#### GTF annotation
raw_annotation_gft_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/gff3_annotations_for_igv/gff3_filtered.gff3"
raw_annotation_gtf <- rtracklayer::import(con = raw_annotation_gft_file)
raw_annotation_gtf <- raw_annotation_gtf[mcols(raw_annotation_gtf)$transcript_id %in% annotation_of_genes_select$ensembl_transcript_id_version]
raw_annotation_gtf_df <- as.data.frame(raw_annotation_gtf)
raw_annotation_gtf_df <- raw_annotation_gtf_df %>% 
  dplyr::select(seqnames, start, end, width,strand, 
                type, ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number) %>% 
  dplyr::filter(transcript_id %in% annotation_of_genes_select$ensembl_transcript_id_version)

remove(ensembl);remove(filts);remove(attrs)
remove(mart);remove(annotation_of_genes);remove(raw_annotation_gtf)

##### Do the plots
# flame_results <- flame_results %>% 
#   dplyr::filter(sample_name == "LRS4_LRS4_barcode02") %>% ### OJO!
#   dplyr::filter(gene_symbol == "CHEK2")

for (s in 1:length(unique(flame_results$sample_name))) {
  #### Data now
  #### Per sample
  sample_now <- unique(flame_results$sample_name)[s]
  data_sample <- flame_results %>% 
    dplyr::filter(sample_name == sample_now)
  
  data_now_final_folder <- paste0(flame_plots_dir,"/",sample_now)
  
  if (!dir.exists(data_now_final_folder)) {
    dir.create(path = data_now_final_folder)
    }
  
  cat(paste0("Working on sample ",s," out of ",length(unique(flame_results$sample_name))," named ",sample_now,"\n"))
  
  #### Per gene
  for (i in 1:length(unique(data_sample$gene_symbol))) {
    
    ## Data of the gene
    gene_now <- unique(data_sample$gene_symbol)[i]
    cat(paste0("Working on gene ",i," out of ",length(unique(data_sample$gene_symbol))," named ",gene_now,"\r"))
    
    data_now <- data_sample %>% 
      dplyr::filter(gene_symbol == gene_now) %>% 
      dplyr::filter(!type %in% c("gene","transcript")) %>% 
      dplyr::mutate(rank = as.numeric(rank),
                    overlap_width = as.numeric(overlap_width),
                    percent_overlap = as.numeric(percent_overlap))
    ## Supporting counts of each of the transcripts
    supp_counts <- data_sample %>% 
      dplyr::filter(type %in% c("transcript"))  %>% 
      dplyr::filter(gene_symbol == gene_now) %>% 
      dplyr::select(source,type,transcript,support_count) %>% 
      dplyr::distinct()
    
    ## Refine the type label and raw_id
    data_now <- data_now %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(
        # type = ifelse(type != "Unknown", gsub(pattern = "", replacement = "", x = raw_id),raw_id)
        type = gsub(pattern = "_[0-9]{1,2}$", replacement = "", x = raw_id)
        ) %>% 
      dplyr::ungroup() %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(raw_id = gsub(pattern = "_", replacement = " ", x = raw_id),
                    raw_id = tools::toTitleCase(raw_id)) %>% 
      dplyr::ungroup()
    
    ## Annotation of the gene from te GTF file
    annotation_of_gene_now <- raw_annotation_gtf_df %>%
      dplyr::filter(gene_name == gene_now) %>%
      dplyr::mutate(type = as.character(type)) %>%
      dplyr::filter(!type %in% c("transcript")) %>%
      dplyr::mutate(sample_name = sample_now,
                    source = "GTF_annotation") %>%
      dplyr::rename(gene_symbol = gene_name,
                    transcript = transcript_id) %>%
      dplyr::mutate(support_count = NA,
                    exon_range = "Whole gene",
                    Parent = NA,
                    overlap_width = NA,
                    percent_overlap = "",
                    raw_id = "GTF_annotation",
                    id = "GTF_annotation",) %>%
      dplyr::rename(rank = exon_number) %>%
      dplyr::select(-c(exon_id))
    
    ## Later on interesting data
    chr_now <- unique(data_now$seqnames) # Chromosome
    
    gene_id_now <- unique(annotation_of_gene_now$gene_id) # Ensembl Gene ID
    gene_id_now_trans <- annotation_of_genes_select %>% # Ensembl Trans ID 
      dplyr::filter(hgnc_symbol == gene_now) %>% 
      dplyr::pull(ensembl_transcript_id_version)
    gene_id_nm_now <- annotation_of_genes_select %>% # NM Gene ID
      dplyr::filter(hgnc_symbol == gene_now) %>% 
      dplyr::pull(transcript_mane_select) 
    
    ## Prepare to bind GTF Annotation and FLAMES data
    data_now <- data_now %>% 
      dplyr::select(-c(diff_start,diff_end,info_from))
    annotation_of_gene_now <- annotation_of_gene_now %>% 
      dplyr::mutate(original_start = start,
                    original_end = end)
    
    all(colnames(annotation_of_gene_now) %in% colnames(data_now))
    annotation_of_gene_now <- annotation_of_gene_now[,
                                                     match(colnames(data_now), colnames(annotation_of_gene_now)),
                                                     drop = FALSE]
    all(colnames(annotation_of_gene_now) %in% colnames(data_now))
    
    #### Generate the data to plot
    data_to_plot <- rbind(annotation_of_gene_now, data_now)
    
    ## Remove not useful columns
    data_to_plot <- data_to_plot %>%
      dplyr::select(-c(sample_name, seqnames, ID, gene_id,Parent))
    
    ### Calculate the introns information (based on original start and end)
    introns <- data_to_plot %>%
      dplyr::filter(type == "exon") %>%
      dplyr::group_by(transcript) %>%
      dplyr::arrange(start, .by_group = TRUE) %>%
      dplyr::mutate(intron_start = lag(original_end) +1, 
                    intron_end = original_start-1) %>%
      dplyr::ungroup() %>% 
      dplyr::filter(!is.na(intron_start)) %>%
      dplyr::relocate(intron_start, .after = original_start) %>%
      dplyr::relocate(intron_end, .after = original_end) %>%
      dplyr::select(-c(start, end)) %>%
      dplyr::rename(start = intron_start, end = intron_end) %>%
      dplyr::mutate(type = "intron")
    
    ## Add intron information
    data_to_plot <- rbind(data_to_plot, introns)
    
    ## Calculate the CDS information 
    data_to_plot_cds <- data_to_plot %>% 
      dplyr::filter(type == "CDS") %>% 
      dplyr::filter(source != "GTF_annotation") %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(cds_rank = as.numeric(x = gsub(pattern = "CDS", replacement = "", x = raw_id))) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(transcript) %>% 
      dplyr::filter(cds_rank == max(cds_rank, na.rm = TRUE) |
                      cds_rank == min(cds_rank, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-cds_rank)
    
    ## Add CDS information
    data_to_plot <- data_to_plot %>%
      dplyr::filter(type != "CDS")
    data_to_plot <- rbind(data_to_plot, data_to_plot_cds)
    
    ### Format of the data to plot information
    data_to_plot <- data_to_plot %>%
      ## ordering of the data
      dplyr::mutate(
        source = factor(source, levels = c("FLAMES", "reference", "GTF_annotation"))) %>%
      dplyr::arrange(source, transcript, start) %>%
      dplyr::mutate(
        ordering = paste0(source, " ", transcript),
        ordering = factor(ordering, levels = unique(ordering))) %>% 
      ## Middle point
      dplyr::rowwise() %>% 
      dplyr::mutate(middle_point = mean(x = c(original_start, original_end))) %>% 
      dplyr::ungroup() %>%
      ## Is it whole
      dplyr::mutate(percent_overlap = ifelse(!is.na(percent_overlap),as.numeric(percent_overlap),percent_overlap)) %>% 
      dplyr::mutate(is_it_whole = ifelse(percent_overlap == 1,"Whole","Non-whole")) %>%
      ## Raw ID formatting
      dplyr::rowwise() %>% 
      dplyr::mutate(raw_id = ifelse(grepl(pattern = "UTR", x = raw_id),
                                    gsub(pattern = " [0-9]{1,3}", replacement = "", x = raw_id),
                                    raw_id)) %>%
      dplyr::mutate(raw_id = ifelse(is_it_whole == "Non-whole",
                                    paste0(raw_id,": ",percent_overlap),
                                    raw_id)) %>%
      dplyr::ungroup()
    
    ### Format the data to plot information for non-whole ids
    data_to_plot_non_complete <- data_to_plot %>% 
      dplyr::filter(source != "GTF_annotation") %>% 
      dplyr::filter(type != "intron") %>%   
      dplyr::filter(is_it_whole == "Non-whole")
    
    ### Format of the supporting counts information (ordering of the data)
    supp_counts <- supp_counts %>%
      dplyr::distinct()
    supp_counts <- supp_counts %>%
      dplyr::mutate(
        source = factor(source, levels = c("FLAMES", "reference", "GTF_annotation"))
      ) %>%
      dplyr::arrange(source, transcript) %>%
      dplyr::mutate(ordering = paste0(source," ",transcript),
                    ordering = factor(ordering, levels = unique(data_to_plot$ordering)))
    
    #### Do the plot
    igv_like_plot <- ggplot(data = data_to_plot) +
      #### Intronic data
      geom_segment(data = . %>%
                     filter(type == "intron"),
                   mapping = aes(x = start, xend = end,
                                 y = ordering, yend = ordering),
                   linewidth = 0.6, colour = "black") +
      #### Exonic data
      ## Segments - original data
      geom_segment(data = . %>%
                     filter(!type %in% c("intron","CDS","Unknown")),
                   mapping = aes(x = original_start, xend = original_end,
                                 y = ordering, yend = ordering),
                   linewidth = 4, colour = "blue") +
      ## Segments - data detected by FLAMES
      geom_segment(data = . %>% 
                     filter(!type %in% c("intron","CDS","Unknown")) %>%
                     filter(!raw_id %in% c("GTF_annotation")) %>%
                     dplyr::filter(is_it_whole == "Non-whole"),
                   mapping = aes(x = start, xend = end, 
                                 y = ordering, yend = ordering), 
                   alpha = 0.5, colour = "green", linewidth = 4) +
      geom_segment(data = . %>% 
                     filter(type %in% c("intron")) %>%
                     filter(!raw_id %in% c("GTF_annotation")) %>%
                     dplyr::filter(is_it_whole == "Non-whole"),
                   mapping = aes(x = start, xend = end, 
                                 y = ordering, yend = ordering), 
                   alpha = 0.5, colour = "pink", linewidth = 5) +
      
      ## Text
      geom_text_repel(data = . %>%
                        filter(!type %in% c("intron","Unknown")) %>%
                        filter(!source %in% c("GTF_annotation")) %>%
                        dplyr::select(ordering,middle_point,raw_id,is_it_whole) %>%
                        filter(!is.na(middle_point)),
                      mapping = aes(x = middle_point, y = ordering,
                                    label = raw_id, colour = is_it_whole),
                      nudge_x = -0.5, nudge_y = -0.5, max.overlaps = Inf) +
      scale_colour_manual(values = c("Whole" = "black",
                                     "Non-whole" = "red")) +
      #### Unknown data
      geom_segment(data = . %>%
                     filter(type == "Unknown"),
                   mapping = aes(x = start, xend = end,
                                 y = ordering, yend = ordering),
                   linewidth = 0.6, colour = "red") +
      #### CDS data
      ## Segments
      geom_segment(data = . %>%
                     filter(type %in% c("CDS")),
                   mapping = aes(x = original_start, xend = original_end,
                                 y = ordering, yend = ordering),
                   linewidth = 2, colour = "purple") +
      #### Text annotations
      geom_label(data = . %>%
                   dplyr::select(ordering, exon_range) %>%
                   distinct(),
                 mapping = aes(x = -Inf, y = ordering,
                               label = paste0("Exon range: ",exon_range)),
                 vjust = -1, hjust = 0) +
      geom_label(data = supp_counts %>%
                   dplyr::select(ordering,support_count) %>%
                   distinct(),
                 mapping = aes(x = Inf, y = ordering,
                               label = paste0("Supp counts: ",support_count)),
                 vjust = -1, hjust = 2) +
      #### Themes
      theme_classic() +
      labs(x = "Position", y = "Transcript", 
           title = paste0(sample_now," - ",gene_now," - ",gene_id_now),
           subtitle = paste0(gene_id_now_trans," used as reference for isoform detection with NM: ",gene_id_nm_now))+
      theme(panel.background = element_rect(fill = "gray"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none")
    
    
    #### Save the plot
    ggsave(filename = paste0(data_now_final_folder,"/",gene_now,"_",sample_now,".png"), plot = igv_like_plot, 
           width = 15, height = 15)
    ggsave(filename = paste0(data_now_final_folder,"/",gene_now,"_",sample_now,".pdf"), plot = igv_like_plot, 
           width = 15, height = 15)
    }
}





