#########################
### Generate manifest ###
#########################

### libraries and WD
setwd("/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/config")
library(dplyr)

fastq_place <- "/imppc/labs/eclab/ijarne/0_Recerca/5_MINION_ENIGMA/concatenated_fastq_files"

### data 
epp <- readr::read_table(file = "./reads_manifest_ENIGMA.txt",
                         col_names = c("sample","cond","batch","fastq"))
epp <- epp %>% 
  dplyr::select(-c(fastq))

epp <- epp %>% 
  dplyr::mutate(fastq = paste0(
    ## base path 
    fastq_place,"/",
    ## folder batch_batch_barcodeXX
    gsub(pattern = "^.*\\.", replacement = "", x = batch),"_",
    gsub(pattern = "^.*\\.", replacement = "", x = batch),"_",
    gsub(pattern = "\\..*$", replacement = "", x = batch),"/",
    ## fastq file batch_batch_barcodeXX.fastq.gz
    gsub(pattern = "^.*\\.", replacement = "", x = batch),"_",
    gsub(pattern = "^.*\\.", replacement = "", x = batch),"_",
    gsub(pattern = "\\..*$", replacement = "", x = batch),".fastq.gz"
    )) %>% 
  dplyr::mutate(exists = file.exists(fastq)) %>% 
  dplyr::select(-c(exists))

readr::write_tsv(x = epp, file = "./reads_manifest_ENIGMA_def.txt", col_names = FALSE)
