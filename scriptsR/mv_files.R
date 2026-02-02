############################
### mv files for MINION ###
############################

### WD and library
setwd("/imppc/labs/eclab/Raw_Data/MINION_IGNASI")
final_place <- "/imppc/labs/eclab/Raw_Data/MINION_IGNASI/AA_final_place"

library(dplyr)

folders <- c("batch_PKD","Prova17","Prova20","Run1_1711")
folders <- c("LRS4","LRS5","LRS6","LRS7","LRS8")

### pick the sample folder
for (folder in 1:length(folders)) {
  
  names_now <- folders[folder]
  folder_now <- paste0("/imppc/labs/eclab/Raw_Data/MINION_IGNASI/",names_now,"/fastq_pass")
  raw_list_of_files <- list.files(path = folder_now,
                                  pattern = ".fastq\\.gz$",
                                  recursive = TRUE, full.names = TRUE)
  cat(paste0("Working on batch ",folder_now," ",folder," out of ",length(folders),"\n"))
  for (i in 1:length(raw_list_of_files)) {
    cat("Working on file ", i, " out of ", length(raw_list_of_files),"\r")
    ## which file (raw_name)
    file_now <- raw_list_of_files[i]
    ## which batch
    batch_name_now <- gsub(pattern = "\\/fastq_pass\\/.*$", replacement = "", x = file_now)
    batch_name_now <- gsub(pattern = ".*/", replacement = "", x = batch_name_now)
    ## which barcode
    barcode_now <- gsub(pattern = ".*\\/(LRS[0-9]{1}_barcode[^/]+\\.*)", replacement = "\\1", x = file_now)
    barcode_now <- gsub(pattern = "(LRS[0-9]{1}_barcode[0-9]{1,3}).*", replacement = "\\1", x = barcode_now)
    
    ## generate the new sample_name
    new_sample_name <- paste0(batch_name_now,"_",barcode_now)
    final_folder <- paste0(final_place,"/",new_sample_name)
    
    ## cmd line1 (generate the folder)
    cmd_1 <- paste0("mkdir -p ",final_folder)
    system(cmd_1)
    
    cmd_2 <- paste0("cp ",file_now," ",final_folder)
    system(cmd_2)
    }
}



