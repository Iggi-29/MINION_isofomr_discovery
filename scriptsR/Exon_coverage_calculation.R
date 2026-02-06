#############################################
### Generate .bed files to check coverage ###
#############################################

#### libraries
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(ggplot2)

#### Constants
### genes to work on
genes_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config/keep_this.txt"
### GTF file
raw_annotation_gtf_file <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run11/ref_files/filtered.gff3"
### list of bams
bams_loc1 <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run11/results/flame"
bams_loc2 <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/results/flame"
### Bed file location
bed_files <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Runs_analysis/bed_files"
### Cov file location
cov_files <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Runs_analysis/coverage"


######### HERE PUT THE SNAKEMAKE PARAMS SO ALL THAT IS ABOVE THIS PART SHOULD BE TAKEN BY SNAKEMAKE
######
###

##### Genes to work on
genes_to_work_on <- readLines(con = genes_to_work_on)

##### BiomaRt work
#### Connect to the service
mart <- biomaRt::useEnsembl(biomart = "ensembl", version = 115)
datasets <- biomaRt::listDatasets(mart = mart)
ensembl <- biomaRt::useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
#### Get the filters and attributes
filts <- biomaRt::listFilters(mart = ensembl)
attrs <- biomaRt::listAttributes(mart = ensembl)
#### Annotations of genes select
annotation_of_genes_select <- biomaRt::getBM(
  mart = ensembl,
  attributes = c("hgnc_symbol","ensembl_transcript_id_version","transcript_mane_select"),
  filters = c("hgnc_symbol"), 
  values = genes_to_work_on)
annotation_of_genes_select <- annotation_of_genes_select %>% 
  dplyr::mutate(is_MANE_select = ifelse(transcript_mane_select == "","no","yes"))

remove(mart);remove(datasets);remove(ensembl);
remove(filts);remove(attrs)

##### GTF importation
raw_annotation_gtf <- rtracklayer::import(con = raw_annotation_gtf_file)
raw_annotation_gtf <- raw_annotation_gtf[mcols(raw_annotation_gtf)$gene_name %in% genes_to_work_on]
raw_annotation_gtf <- as.data.frame(raw_annotation_gtf) %>% 
  dplyr::select(seqnames,start,end,width,strand,
                type,ID,
                gene_id,gene_name,
                transcript_id,
                exon_id,exon_number)

raw_annotation_gtf <- raw_annotation_gtf %>% 
  dplyr::mutate(type = as.character(type)) %>% 
  dplyr::filter(type %in% "exon") %>% 
  dplyr::filter(transcript_id %in% annotation_of_genes_select$ensembl_transcript_id_version[annotation_of_genes_select$is_MANE_select == "yes"])


for (i in 1:length(unique(raw_annotation_gtf$gene_name))){
 gene_now <- unique(raw_annotation_gtf$gene_name)[i] 
 annotation_now <- raw_annotation_gtf %>% 
   dplyr::filter(gene_name == gene_now)
 annotation_now <- annotation_now %>% 
   dplyr::select(seqnames, start, end, gene_name,exon_number,strand)
 annotation_now <- annotation_now %>% 
   dplyr::mutate(id = paste0(gene_name,"_","exon",exon_number)) %>% 
   dplyr::select(-c(gene_name,exon_number)) %>% 
   dplyr::relocate(id, .after = end) %>% 
   dplyr::mutate(dunno = ".") %>% 
   dplyr::select(seqnames,start,end,id,dunno,strand)
 
 fileName_final <- paste0(x = c(bed_files,"/",gene_now,".bed"),collapse = "")
 readr::write_tsv(x = annotation_now, col_names = FALSE, file = fileName_final)
 }


# barcode04_sorted.bam
### coverage of bam fles
## bams
bam_files1 <- list.files(path = bams_loc1, pattern = "_sorted\\.bam$", full.names = TRUE, recursive = FALSE)
bam_files2 <- list.files(path = bams_loc2, pattern = "_sorted\\.bam$", full.names = TRUE, recursive = FALSE)
bam_files <- c(bam_files1, bam_files2)
remove(bam_files2);remove(bam_files1)

## bed
bed_files <- list.files(path = bed_files, pattern = "\\.bed$", full.names = TRUE, recursive = FALSE)
bed_files <- bed_files[grep(pattern = "NF2|NF1", x = bed_files)]

for (i in 1:length(bam_files)) {
  ### data_now
  bam_file_now <- bam_files[i]
  sample_name_now <- gsub(x = basename(bam_file_now), pattern = "_sorted\\.bam$", replacement = "")
  
  for (b in 1:length(bed_files)) {
    ## gene now
    bed_file_now <- bed_files[b]
    gene_now <- gsub(pattern = "\\.bed$", replacement = "", x = basename(bed_file_now))
    ## prepare the commands
    basic_bed_tools_mean <- paste0(x = c("bedtools coverage -a ",bed_file_now," -b ",bam_file_now," -mean"), collapse = "")
    basic_bed_tools_hist <- paste0(x = c("bedtools coverage -a ",bed_file_now," -b ",bam_file_now," -hist"), collapse = "")
    
    cmd_mean <- paste0(x = c(basic_bed_tools_mean," > ",cov_files,"/",sample_name_now,"_",gene_now,"_depth.txt"), collapse = "")
    cmd_hist <- paste0(x = c(basic_bed_tools_hist," > ",cov_files,"/",sample_name_now,"_",gene_now,"_coverage.txt"), collapse = "")
  
    ## execute the commands
    system(cmd_mean)
    # system(cmd_hist)
    }
  }


remove(annotation_now);remove(annotation_of_genes_select);remove(i);remove(b);
remove(raw_annotation_gtf);remove(bam_file_now);remove(bams_loc1);remove(bams_loc2)
remove(genes_to_work_on);remove(bam_files);remove(basic_bed_tools_hist);remove(basic_bed_tools_mean)
remove(bed_file_now);remove(bed_files);remove(fileName_final)
remove(gene_now);remove(cmd_mean);remove(cmd_hist)
remove(raw_annotation_gtf_file);remove(sample_name_now)

#### Cov file importation
cov_files_place <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Runs_analysis/coverage"

cov_files <- list.files(path = cov_files_place, pattern = "_depth\\.txt$", full.names = TRUE, recursive = FALSE)
cov_files_list <- lapply(X = cov_files, FUN = function(x) {
  readr::read_tsv(file = x, col_names = c("seqnames","start","end","region","dunno","strand","mean_depth")) %>% 
    dplyr::mutate(sample_name = gsub(x = basename(x), pattern = "_depth\\.txt$", replacement = ""))
})
cov_files_list <- do.call(rbind, cov_files_list)
cov_files_list <- cov_files_list %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(gene = unlist(strsplit(x = region, split = "_"))[[1]][1],
                region = unlist(strsplit(x = region, split = "_"))[[2]][1]) %>%
  dplyr::ungroup()


for (i in 1:length(unique(x = cov_files_list$gene))) {
  ### data now
  gene_now <- unique(cov_files_list$gene)[i]
  data_now <- cov_files_list %>% 
    dplyr::filter(gene == gene_now)
  
  data_now <- data_now %>% 
    dplyr::arrange((start)) %>% 
    dplyr::mutate(region = factor(region, levels = unique(region))) %>% 
    dplyr::ungroup()
  
  plot_of_cov <- ggplot(data = data_now) +
    geom_col(mapping = aes(x = region, y = mean_depth))+
    labs(x = "", y = "", title = paste0("Mean coverage of exons for ",gene_now))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(. ~ sample_name, scales = "free_y")
  
  ggsave(filename = paste0(x = c(cov_files_place,"/plots/",gene_now,".png"), collapse = ""), plot = plot_of_cov, width = 20, height = 20)
  ggsave(filename = paste0(x = c(cov_files_place,"/plots/",gene_now,".pdf"), collapse = ""), plot = plot_of_cov, width = 20, height = 20)
  
}









  
  