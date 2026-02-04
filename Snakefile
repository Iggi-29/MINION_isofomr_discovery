#########################################################
### File to organize rules for MINION data processing ###
#########################################################

#### libraries
import pandas as pd
from shutil  import copyfile
import os, glob, re, sys

#### Configuration of the works to do
### config file
configfile: "./config/config.yaml"

### variables definition
## envs
AGAT_anaconda = config["anaconda_envs"]["AGAT_anaconda"]
sequinR_anaconda = config["anaconda_envs"]["sequinR_anaconda"]
fastpqc_env = config["anaconda_envs"]["fastpqc_env"]
bamslam_env = config["anaconda_envs"]["bamslam_env"]
flair_env = config["anaconda_envs"]["flair_env"]
flair_env2 = config["anaconda_envs"]["flair_env2"]
flame_env = config["anaconda_envs"]["flame_env"]
readr_env = config["anaconda_envs"]["readr_env"]
biomart_env = config["anaconda_envs"]["biomart_env"]

## images
minimap2_img = config["singularity_imgs"]["minimap2_img"]
samtools_img = config["singularity_imgs"]["samtools_img"]
fastqc_image = config["singularity_imgs"]["fastqc_image"]
multiqc_image = config["singularity_imgs"]["multiqc_image"]

## manifest
manifest_file = config["manifest_file"]
manifest_file_lrs4 = config["manifest_file_lrs4"]
manifest_file_lrs5 = config["manifest_file_lrs5"]
manifest_file_lrs6 = config["manifest_file_lrs6"]
manifest_file_lrs7 = config["manifest_file_lrs7"]
manifest_file_lrs8 = config["manifest_file_lrs8"]

## samples
samples = config["samples"]
samples_lrs4 = config["samples_lrs4"]
samples_lrs5 = config["samples_lrs5"]
samples_lrs6 = config["samples_lrs6"]
samples_lrs7 = config["samples_lrs7"]
samples_lrs8 = config["samples_lrs8"]
manifest_file_lrs1 = config["manifest_file_lrs1"]

## genes
genes_to_work_on = config["genes_to_work_on"]
genes_to_filter = config["genes_to_filter"]
genes_to_work_on_ensembl = config["genes_to_work_on_ensembl"]
## ref files
genome = config["genome"]
transcriptome = config["transcriptome"]
gtf_annotation = config["gtf_annotation"]
gff3_annotation = config["gff3_annotation"]

## paths
base_path = config["base_path"]
input_files_folder = config["input_files_folder"]
snake_files = config["snake_files"]
scriptsR = config["scriptsR"]
flame_script = config["flame_script"]
flame_config = config["flame_config"] 

## color
colors = config["colors"]

#### Include the rules of the analysis
include: "Prepare_data.smk"
include: "File_concatenation.smk"
include: "QC_report_generation.smk"
include: "Transcriptome_Alignment.smk"
include: "flair.smk"
include: "flame.smk"

### rules
rule all:
    input:
        rules.trigger_file_concatenation.input,
        rules.trigger_prepare_data.input,
        rules.trigger_QC.input,
        rules.trigger_transcriptome_align.input,
        rules.trigger_flame.input,
        # rules.trigger_flair.input

