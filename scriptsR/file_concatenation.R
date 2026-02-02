###########################
### fastQ concatenation ###
###########################

library(dplyr)

### Snakemake
final_gff3 <- snakemake@output[["gff3_filtered_MANE"]]
out_dir <- snakemake@params[["out_dir"]]
input_file_place <- snakemake@params[["input_file_place"]]
samples <- snakemake@params[["samples"]]

### trial
# out_dir <- "/imppc/labs/eclab/ijarne/0_Recerca/6_Run2/concatenated_fastq_files"
# input_file_place <- "/imppc/labs/eclab/Raw_Data/MINION_IGNASI/Run2_14012026/no_sample_id/20260114_1249_MC-114182_FBE54097_668d5d39/fastq_pass"
# samples <- c("barcode03","barcode13","barcode14")

### Do the concatenation
for (i in 1:length(samples)) {
  cat(paste0(i,"\n"))
  ### data now
  sample_now <- samples[i]
  # input_now <- paste0(
  #   x = c(input_file_place,"/",sample_now,"/"),
  #   collapse = "")
  input_now <- paste0(
    x = c(input_file_place,"/"),
    collapse = "")


  folder_out_now <- paste0(
    x = c(out_dir,"/",sample_now), 
    collapse = ""
    )
  final_out_name <- paste0(
    x = c(out_dir,"/",sample_now,"/",sample_now,".fastq.gz"), 
    collapse = "")
  
  list_of_files_to_concat <- list.files(path = input_now, pattern = "\\.fastq\\.gz", full.names = TRUE)
  cat(list_of_files_to_concat)

  ### do the cmd work
  cmd_one <- paste0(x = c("mkdir -p ",folder_out_now), collapse = "")
  system(cmd_one)
  
  cmd_two <- paste0(
    x = c("cat ",paste0(x = list_of_files_to_concat, collapse = " "),
          " > ",shQuote(final_out_name)), collapse = "")
  cat(cmd_two)
  system(cmd_two)
  
  }













