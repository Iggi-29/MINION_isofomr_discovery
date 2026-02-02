#################################
### FastQC file concatenation ###
#################################

### fastQ file concatenation
rule fastq_concatenation:
    input:
        "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/start_the_pipeline.txt"
    output:
        concatenated_fastq_files = base_path + "/concatenated_fastq_files/{sample}/{sample}.fastq.gz"
    params:
        out_dir = base_path + "/concatenated_fastq_files/",
        input_file_place = input_files_folder + "/{sample}",
        samples = "{sample}",
        scriptsR = scriptsR
    threads:
        2
    log:
        terminal_log = snake_files + "/log/{sample}_fastq_concatenation.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_fastq_concatenation.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_fastq_concatenation.bmk"
    script:"{params.scriptsR}/file_concatenation.R"

rule trigger_file_concatenation:
    input:
        expand(base_path + "/concatenated_fastq_files/{sample}/{sample}.fastq.gz",
        base_path = base_path, sample = samples),

