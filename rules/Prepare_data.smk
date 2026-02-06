###############################
### Prepare the data itself ###
###############################

# ### fastQ file concatenation
# rule fastq_concatenation:
#     input:
#         "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/MINION/start_the_pipeline.txt"
#     output:
#         concatenated_fastq_files = base_path + "/concatenated_fastq_files/{sample}/{sample}.fastq.gz"
#     params:
#         out_dir = base_path + "/concatenated_fastq_files/{sample}/",
#         input_file_place = input_files_folder + "/{sample}"
#     threads:
#         2
#     log:
#         terminal_log = snake_files + "/log/{sample}_fastq_concatenation.log", 
#         snakemake_log = snake_files + "/snakemake_log/{sample}_fastq_concatenation.log"
#     benchmark: 
#         snake_files + "/benchmark/{sample}_fastq_concatenation.bmk"
#     shell:"""
#         (
#             mkdir -p "$(dirname {output.concatenated_fastq_files})"

#             files=({params.input_file_place}/*.fastq.gz)

#             if [ -e "${{files[0]}}" ]; then
#                 cat {params.input_file_place}/*.fastq.gz > "{output.concatenated_fastq_files}"
#                 echo "Created {output.concatenated_fastq_files}"
#             else
#                 echo "Warning: no FASTQ files found in {params.input_file_place}"
#                 touch "{output.concatenated_fastq_files}"
#             fi
#         ) > "{log.terminal_log}" 2>> "{log.snakemake_log}"
#         """


### GFF3 file filtration
rule gff3_file_filtration:
    input:
        genes_to_filter = genes_to_filter,
        gff3_annotation = gff3_annotation
    output:
        gff3_filtered = base_path + "/ref_files/filtered.gff3"
    conda:
        AGAT_anaconda
    log:
        terminal_log = snake_files + "/log/gff3_file_filtration.log", 
        snakemake_log = snake_files + "/snakemake_log/gff3_file_filtration.log"
    benchmark: 
        snake_files + "/benchmark/gff3_file_filtration.bmk"
    shell:"""
        agat_sp_filter_feature_from_keep_list.pl -f {input.gff3_annotation} \
        -a gene_name -kl {input.genes_to_filter} \
        -o {output.gff3_filtered} > {log.terminal_log} 2>> {log.snakemake_log}
        """

### GFF3 file filtration only MANE
rule gff3_file_filtration_MANE:
    input:
        gff3_filtered = base_path + "/ref_files/filtered.gff3"
    output:
        gff3_filtered_MANE = base_path + "/ref_files/filtered_MANE.gff3"
    conda:
        readr_env
    params:
        scriptsR = scriptsR
    log:
        terminal_log = snake_files + "/log/gff3_file_filtration_MANE.log", 
        snakemake_log = snake_files + "/snakemake_log/gff3_file_filtration_MANE.log"
    benchmark: 
        snake_files + "/benchmark/gff3_file_filtration.bmk"
    script:"{params.scriptsR}/gff3_filtration.R"


### Transcriptome filtering
rule transcriptome_filtration:
    input:
        transcriptome = transcriptome
    params:
        scriptsR = scriptsR,
        genes_to_work_on = "{genes_to_work_on}"
    output:
        transcriptome_filtered = base_path + "/ref_files/transciptome_filtered/{genes_to_work_on}.fa"
    conda:
        sequinR_anaconda
    log:
        terminal_log = snake_files + "/log/{genes_to_work_on}_transcriptome_filtration.log", 
        snakemake_log = snake_files + "/snakemake_log/{genes_to_work_on}_transcriptome_filtration.log"
    benchmark: 
        snake_files + "/benchmark/{genes_to_work_on}_transcriptome_filtration.bmk"
    script:"{params.scriptsR}/filter_Transcriptome.R"

rule trigger_prepare_data:
    input:
        # expand(base_path + "/concatenated_fastq_files/{sample}/{sample}.fastq.gz",
        # base_path = base_path, sample = samples),
        base_path + "/ref_files/filtered.gff3",
        base_path + "/ref_files/filtered_MANE.gff3",
        expand(base_path + "/ref_files/transciptome_filtered/{genes_to_work_on}.fa",
        base_path = base_path, genes_to_work_on = genes_to_work_on)

