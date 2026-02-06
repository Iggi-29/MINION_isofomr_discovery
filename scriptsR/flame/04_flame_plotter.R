#######################################
### FLAMES long reads data plotting ###
#######################################

#### libraries
library(ggplot2)
library(ggnewscale)
library(ggrepel)

library(openxlsx)
library(readxl)

library(dplyr)

#### Constants
### flame plots dir
flame_plots_dir <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flame/final_select/plots"
### flame results exons
flame_data <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flame/final_select/matched_exon_data.xlsx"
flame_data_gene_and_trans <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flame/final_select/matched_gene_and_trans_data.xlsx"
### trans annotation
annotaton_of_trans <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flame/final_select/annotation_of_transcripts.xlsx"
### Events
event_df <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA2/results/flame/final_select/events.xlsx"

######### HERE PUT THE SNAKEMAKE PARAMS SO ALL THAT IS ABOVE THIS PART SHOULD BE TAKEN BY SNAKEMAKE
######
###

#### flame plots dir
if (!dir.exists(paths = flame_plots_dir)) {
  dir.create(path = flame_plots_dir)
  cat(paste0("Directory to save FLAME plots has been created","\n"))
} else {
  cat(paste0("Directory to save FLAME plots exists, no further action is taken.","\n"))
}

#### Import the data
### flame data
sample_names <- readxl::excel_sheets(path = flame_data)
flame_results <- lapply(X = sample_names, FUN = function(x) {
  openxlsx::read.xlsx(xlsxFile = flame_data, sheet = x)
})
names(flame_results) <- sample_names

### flame data exons and trans
sample_names <- readxl::excel_sheets(path = flame_data_gene_and_trans)
flame_results_gene_and_trans <- lapply(X = sample_names, FUN = function(x) {
  openxlsx::read.xlsx(xlsxFile = flame_data_gene_and_trans, sheet = x)
})
names(flame_results_gene_and_trans) <- sample_names

### annotation data
annotaton_of_trans <- openxlsx::read.xlsx(xlsxFile = annotaton_of_trans)

### event data
event_df <- openxlsx::read.xlsx(xlsxFile = event_df)
all(event_df$sample_name %in% c(names(flame_results), names(flame_results_gene_and_trans)))

for (i in 1:length((event_df$sample_name))) {
  ### data now
  ## event data
  sample_now <- unique(event_df$sample_name)[i]
  event_df_now <- event_df %>% 
    dplyr::filter(sample_name == sample_now)
  event_now_gene <- event_df_now$gene_symbol
  ## flame data
  # exons
  flame_results_now <- flame_results[[sample_now]] %>% 
    dplyr::filter(gene_symbol == event_now_gene) %>% 
    dplyr::mutate(source = ifelse(grepl(x = transcript_id, pattern = "ENST"),"reference","FLAME"))
  # gene and transcrpt data
  flame_results_gene_and_trans_now <- flame_results_gene_and_trans[[sample_now]] %>% 
    dplyr::filter(is.na(Parent) | raw_id %in% (flame_results_now$Parent)) %>% 
    dplyr::filter((!is.na(ensembl_gene_id) & ensembl_gene_id %in% flame_results_now$ensembl_gene_id) | raw_id %in% (flame_results_now$Parent)) %>% 
    dplyr::select(strand,ensembl_gene_id,support_count,info_from,Parent,raw_id) %>% 
    dplyr::mutate(gene_symbol = event_now_gene) %>% 
    dplyr::select(gene_symbol, ensembl_gene_id, strand, info_from, Parent, support_count,raw_id) %>% 
    dplyr::mutate(id = paste0(sample_now,"_",gene_symbol,"_",raw_id))
  
  ## annotation data
  annotation_of_trans_now <- annotaton_of_trans %>%
    dplyr::filter(gene_name == event_now_gene) %>%
    dplyr::filter(transcript_id %in% gsub(x = flame_results_now$reff_transcript, pattern = "transcript\\:", replacement = "") |
                    is_MANE_select == "yes")
  
  annotation_of_trans_now <- annotation_of_trans_now %>% 
    dplyr::rename(original_start = start,
                  original_end = end,
                  raw_id = ID,
                  ensembl_gene_id = gene_id,
                  gene_symbol = gene_name) %>% 
    dplyr::filter(type %in% c("exon","intron")) %>% 
    dplyr::mutate(sample_name = sample_now,
                  start = original_start,
                  end = original_end,
                  diff_start = start,
                  diff_end = end,
                  source = "GTF_annotation",
                  reff_transcript = "",ensembl_coincidence = "",exon_range = "",Parent = "", 
                  rank = NA, overlap_width = NA, percent_overlap = NA,
                  ensembl_coincidence = "",
                  info_from = type, support_count = NA, id = NA, exon_number = NA, is_MANE_select = "") %>% 
    dplyr::rename(ensembl_coincidence_id = ensembl_coincidence) %>% 
    dplyr::mutate(id = paste0(sample_name,"_",gene_symbol,"_transcript:",transcript_id))
  
  annotation_of_trans_now <- annotation_of_trans_now[,grep(pattern = paste0(x = colnames(flame_results_now), collapse = "|"), 
                                                           x = colnames(annotation_of_trans_now)), drop = FALSE]
  
  setdiff(x = colnames(annotation_of_trans_now), y = colnames(flame_results_now))
  setdiff(x = colnames(flame_results_now), y = colnames(annotation_of_trans_now))
  
  ### Generate data to plot
  data_to_plot <- rbind(annotation_of_trans_now, flame_results_now)
  # ## order the data
  # data_to_plot <- data_to_plot %>% 
  #   dplyr::mutate(source = factor(source, levels = c("FLAME","reference","Raw_annotation"))) %>% 
  #   dplyr::mutate(ordering = paste0(as.character(source),"_",transcript_id),
  #                 ordering = factor(ordering, levels = unique(ordering)))

  #### Prepare the data to plot
  ## order the data
  data_to_plot <- data_to_plot %>% 
    dplyr::mutate(
      source = factor(source, levels = c("GTF_annotation", "reference", "FLAME"))) %>% 
    dplyr::arrange(source, transcript_id) %>%
    dplyr::mutate(
      ordering = paste(source, transcript_id, sep = "_"),
      ordering = factor(ordering, levels = rev(unique(ordering)))
    ) %>%   
    ## middle point
    dplyr::rowwise() %>% 
    dplyr::mutate(middle_point = mean(x = c(original_start, original_end))) %>% 
    dplyr::ungroup() %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(middle_point = ifelse(type == "extent", mean(x = c(start, end)), middle_point)) %>% 
    dplyr::ungroup() %>% 
    ## is it whole
    # Annotation of coding things
    dplyr::mutate(is_it_whole =
                    ifelse((!type %in% c("extent","intron")) & (percent_overlap != 1) & (as.character(source) != "GTF_annotation"),
                           "Non-whole",
                           ifelse((!type %in% c("extent","intron")) & (percent_overlap == 1) & (as.character(source) != "GTF_annotation"),
                                  "Whole","EPP"))) %>%
    ## Annotation of non-coding things
    dplyr::mutate(is_it_whole = 
                    ifelse((type %in% c("extent","intron")) & (percent_overlap != 0 & as.character(source) != "GTF_annotation"),
                           "Non-whole",
                           ifelse((type %in% c("extent","intron")) & (percent_overlap == 0 & as.character(source) != "GTF_annotation"),
                                  "Whole",is_it_whole))) %>% 
    ## Raw ID format
    dplyr::ungroup() %>% 
    ## NO present work
    dplyr::rowwise() %>% 
    dplyr::mutate(
      start = ifelse((exon_id == "NO PRESENT") & (type == "exon"),
                     original_start, start),
      end = ifelse((exon_id == "NO PRESENT") & (type == "exon"),
                   original_end, end)) %>%
    dplyr::ungroup() %>% 
    ## ID to annotate
    dplyr::rowwise() %>% 
    dplyr::mutate(abs_diff = {
      width_original = original_start - original_end
      width_found = abs(start - end)
    }) %>% 
    dplyr::mutate(annotation_id = paste0(raw_id,"_",percent_overlap,"_",abs_diff)) %>% 
    dplyr::ungroup()
  
  #### Prepare the data to plot
  data_to_plot_trans <- flame_results_gene_and_trans_now %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(source = ifelse((info_from == "gene") & (grepl(pattern = "ENSG", x = raw_id)),"GTF_annotation",
                                  ifelse((info_from != "gene") & (grepl(pattern = "ENST", x = raw_id)),"reference",
                                         ifelse((info_from != "gene") & (grepl(pattern = "ENSG", x = raw_id)),"FLAME",NA)))) %>%
    dplyr::mutate(raw_id = gsub(pattern = "gene\\:", x = raw_id, replacement = "")) %>%
    dplyr::mutate(ordering = paste0(source,"_",raw_id)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(info_from != "gene")
  
  #### The plot
  igv_plot <- ggplot(data = data_to_plot) +
    ### Intronic data
    geom_segment(data = . %>% 
                   filter(type == "intron" & as.character(source) != "GTF_annotation"), 
                 mapping = aes(x = original_start, xend = original_end,
                               y = ordering, yend = ordering)) +
    geom_segment(data = . %>% 
                   filter(type == "intron" & as.character(source) != "GTF_annotation") %>% 
                   filter(is_it_whole == "Non-whole"),
                 mapping = aes(x = start, xend = end,
                               y = ordering, yend = ordering), colour = "firebrick1", linewidth = 5, alpha = 0.75) +
    geom_text_repel(data = . %>% 
                      filter(type == "intron" & as.character(source) != "GTF_annotation") %>% 
                      filter(is_it_whole == "Non-whole"), 
                    mapping = aes(x = middle_point, y = ordering,
                                  label = raw_id), colour = "firebrick1") +
    ### Extent data
    geom_segment(data = . %>% 
                   filter(type == "extent" & as.character(source) != "GTF_annotation"),
                 mapping = aes(x = start, xend = end,
                               y = ordering, yend = ordering), colour = "darkgreen", linewidth = 5, alpha = 0.75) +
    geom_text_repel(data = . %>% 
                      filter(type == "extent" & as.character(source) != "GTF_annotation"), 
                    mapping = aes(x = middle_point, y = ordering,
                                  label = annotation_id), colour = "darkgreen") +
    ### Exonic data
    geom_segment(data = . %>% 
                   filter(type == "exon" & as.character(source) != "GTF_annotation"), 
                 mapping = aes(x = original_start, xend = original_end,
                               y = ordering, yend = ordering), linewidth = 5, colour = "darkblue") +
    geom_segment(data = . %>%
                   filter(type == "exon" & as.character(source) != "GTF_annotation") %>%
                   filter(is_it_whole == "Non-whole"),
                 mapping = aes(x = start, xend = end,
                               y = ordering, yend = ordering), colour = "purple", linewidth = 5, alpha = 1) +
    geom_text_repel(data = . %>% 
                      filter(type == "exon" & as.character(source) != "GTF_annotation") %>% 
                      filter(is_it_whole == "Non-whole"), 
                    mapping = aes(x = middle_point, y = ordering,
                                  label = annotation_id), colour = "purple") +
    #### GTF data
    geom_segment(data = . %>%
                   filter(type == "intron" & as.character(source) == "GTF_annotation"),
                 mapping = aes(x = original_start, xend = original_end,
                               y = ordering, yend = ordering)) +
    geom_segment(data = . %>%
                   filter(type == "exon" & as.character(source) == "GTF_annotation"),
                 mapping = aes(x = original_start, xend = original_end,
                               y = ordering, yend = ordering), linewidth = 5, colour = "darkblue") +
    #### Text annotation
    geom_label(data = data_to_plot_trans,
               mapping = aes(x = Inf, y = ordering, 
                             label = paste0("Supp counts: ",support_count),
                             vjust = -1, hjust = 2)) +
    #### Theme adjustments
    theme_classic() +
    labs(x = "Position", y = "Transcript",
         title = paste0(sample_now," - ",event_now_gene),
         subtitle = paste0("Total reads to the gene: ",unique(flame_results_gene_and_trans_now$support_count[flame_results_gene_and_trans_now$info_from == "gene"]))) +
    theme(panel.background = element_rect(fill = "gray"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  ggsave(plot = igv_plot, filename = paste0(flame_plots_dir,"/",sample_now,"_",event_now_gene,".pdf"), width = 15, height = 15)
  
  }





