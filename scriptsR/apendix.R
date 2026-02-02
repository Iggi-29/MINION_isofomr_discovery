flame_results_df <- flame_results
ranges_to_check <- ranges_select
ranges_corrector <- flame_results_percent_to_match_keep
df_annot <- raw_annotation_gtf_df_select

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
    
    matching_df <- matching_df %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(valid_col = ifelse((grepl(pattern = "UTR|codon|CDS", x = subject_id)) & 
                                         (grepl(pattern = "exon", x = refference_id)),
                                       "NO","YES")) %>% 
      # dplyr::mutate(valid_col = ifelse((grepl(pattern = "UTR|codon|exon", x = subject_id)) & 
      #                                    (grepl(pattern = "CDS", x = refference_id)),
      #                                  "NO",valid_col)) %>%
      dplyr::filter(valid_col == "YES") %>% 
      dplyr::select(-c(valid_col))

# }}
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
    
    ### Remove NO present introns
    matching_df <- matching_df %>% 
      dplyr::rowwise() %>% 
      dplyr::filter((subject_id != "NO present") & (!grepl(pattern = "intron", x = refference_id))) %>% 
      dplyr::ungroup()
    
    ### finish the data.frames
    df_comparisons <- rbind(df_comparisons,matching_df)

  }
}

df_comparisons <- df_comparisons %>%
  ## ID ID match for GTF file based info
  dplyr::rowwise() %>%
  dplyr::mutate(id_id = paste0(refference_id,"_",id)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent_overlap = ifelse(id_id %in% ranges_corrector$id_id,
                                         ranges_corrector$percent_overlap[match(id_id, ranges_corrector$id_id)],
                                         percent_overlap)) %>%
  dplyr::ungroup()




# df_comparisons <- df_comparisons %>% 
#   dplyr::group_by(subject_id,id) %>%
#   dplyr::mutate(subject_id_n = n()) %>% 
#   dplyr::ungroup() %>%
#   dplyr::group_by(refference_id,id) %>% 
#   dplyr::mutate(refference_id_n = n()) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::filter(subject_id != "NO present") %>% 
#   dplyr::filter(subject_id_n > 1 | refference_id_n > 1)

