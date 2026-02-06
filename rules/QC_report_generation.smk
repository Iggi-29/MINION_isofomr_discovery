#################################
###          QC report        ###
#################################

### Fastp
rule fastp:
    input:
        fastq_concat = base_path + "/concatenated_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        json_out = base_path + "/flitered_fastq_files/{sample}/{sample}_report.json",
        html_out = base_path + "/flitered_fastq_files/{sample}/{sample}_report.html",
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    params:
        length_required = 500,
        q = 10,
        cut_front_mean_quality = 8,
        cut_tail_mean_quality = 8
    conda:
        fastpqc_env
    threads:
        2
    log:
        terminal_log = snake_files + "/log/{sample}_fastp.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_fastp.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_fastp.bmk"
    shell:"""
        fastplong --length_required {params.length_required} \
        --thread {threads} \
        --cut_front --cut_tail --cut_front_mean_quality {params.cut_front_mean_quality} \
        --cut_tail_mean_quality {params.cut_tail_mean_quality} \
        --in {input.fastq_concat} --out {output.fastq_filtered} \
        --html {output.html_out} --json {output.json_out} > "{log.terminal_log}" 2>> "{log.snakemake_log}"
        """

### BAMs are like:
## BAM whole <- base_path + "/results/alignments/transcriptome_whole/{sample}.bam"
## BAM per gene <- base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.bam"

### BAMSlam whole transcriptome
rule bam_slam_whole:
    input:
        bam_transcriptome = base_path + "/results/alignments/transcriptome_whole/{sample}.bam"
    output:
        bam_slam_done = base_path + "/results/bamslam_whole_transcriptome/{sample}/bam_slam_done.txt"
    params:
        scriptsR = scriptsR,
        bam_slam_out_name = base_path + "/results/bamslam_whole_transcriptome/{sample}/"
    threads:
        2
    conda:
        bamslam_env
    log:
        base_path + "/Run/snakemake_log/{sample}_bam_slam.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}_bam_slam.bmk"
    shell:"""
    Rscript {params.scriptsR}/BamSlam.R cdna {input.bam_transcriptome} {params.bam_slam_out_name} {output.bam_slam_done}
    """

### BAMSlam per gene
rule bam_slam_per_gene:
    input:
        bam_transcriptome = base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.bam"
    output:
        bam_slam_done = base_path + "/results/bam_slam_per_gene/{sample}/{genes_to_work_on}/{genes_to_work_on}_bam_slam_done.txt"
    params:
        scriptsR = scriptsR,
        bam_slam_out_name = base_path + "/results/bam_slam_per_gene/{sample}/{genes_to_work_on}/{genes_to_work_on}"
    threads:
        2
    conda:
        bamslam_env
    log:
        base_path + "/Run/snakemake_log/{sample}_{genes_to_work_on}_bam_slam.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}_{genes_to_work_on}_bam_slam.bmk"
    shell:"""
    Rscript {params.scriptsR}/BamSlam.R cdna {input.bam_transcriptome} {params.bam_slam_out_name} {output.bam_slam_done}
    """

### FastQC
rule fastQC:
    input:
        reads = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        out_dir = directory(base_path + "/fastQC/{sample}/")
    container:
        fastqc_image
    threads:
        2
    benchmark:
        base_path + "/Run/benchmarks/{sample}_FastQC.bmk"
    log:
        base_path + "/Run/snakemake_log/{sample}_FastQC.log"
    shell:"""
        mkdir {output.out_dir}
        fastqc -t {threads} {input.reads} -o {output.out_dir}  2>> {log}
        """

### MultiQC preparation
# rule multi_QC_input:
#     input:
#         out_dir = expand(base_path + "/fastQC/{sample}/",
#         sample = samples)
#     output:
#         report_list = base_path + "/multiQC/list_of_files.txt",
#     params:
#         base_path = base_path,
#         sub_path = "*/fastQC/*_fastqc.html"
#     threads:
#         2
#     benchmark:
#         base_path + "/Run/benchmarks/MultiQC_input.bmk"
#     log:
#         base_path + "/Run/snakemake_log/MultiQC_input.log"
#     shell:"""
#         dirname $(find {params.base_path} -type f -path '{params.sub_path}') | sort -u > {output.report_list}
#         """

# ### MultiQC
# rule multiQC:
#     input:
#         report_list = base_path + "/multiQC/list_of_files.txt"
#     output:
#         html_report = base_path + "/multiQC/MultiQC_report.html"
#     params:
#         out_dir = base_path + "/multiQC/",
#         html_report = base_path + "/multiQC/MultiQC_report.html"
#     container:
#        multiqc_image
#     threads:
#         2
#     log:
#         base_path + "/Run/benchmarks/MultiQC_input.log"
#     benchmark:
#         base_path + "/Run/benchmarks/MultiQC_input.bmk"
#     shell:"""
#         multiqc --file-list {input.report_list} --filename {params.html_report}
#         """

rule trigger_QC:
    input:
        expand(base_path + "/concatenated_fastq_files/{sample}/{sample}.fastq.gz",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/bamslam_whole_transcriptome/{sample}/bam_slam_done.txt",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/bam_slam_per_gene/{sample}/{genes_to_work_on}/{genes_to_work_on}_bam_slam_done.txt",
        base_path = base_path, sample = samples, genes_to_work_on = genes_to_work_on),
        # html_report = base_path + "/multiQC/MultiQC_report.html"
