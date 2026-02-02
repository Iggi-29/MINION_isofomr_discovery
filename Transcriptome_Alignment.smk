#########################################################
### Alignment against the transcriptome per each gene ###
#########################################################

### Transcriptome alignment per gene
rule transcriptome_alignment:
    input:
        transcriptome_filtered = base_path + "/ref_files/transciptome_filtered/{genes_to_work_on}.fa",
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        sam_transcriptome = temp(base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.sam")
        # sam_transcriptome = base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.sam"
    container:
        minimap2_img
    log:
        terminal_log = snake_files + "/log/{genes_to_work_on}_{sample}_transcriptome_alignment.log", 
        snakemake_log = snake_files + "/snakemake_log/{genes_to_work_on}_{sample}_transcriptome_alignment.log"
    benchmark: 
        snake_files + "/benchmark/{genes_to_work_on}_{sample}_transcriptome_alignment.bmk"
    shell:"""
        minimap2 -ax map-ont --sam-hit-only \
        {input.transcriptome_filtered} \
        {input.fastq_filtered} > {output.sam_transcriptome} 2>> {log.snakemake_log}
        """

rule sam_to_bam_trans:
    input:
        sam_transcriptome = base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.sam"
    output:
        bam_transcriptome = base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.bam"
    container:
        samtools_img
    log:
        terminal_log = snake_files + "/log/{genes_to_work_on}_{sample}_samtools.log", 
        snakemake_log = snake_files + "/snakemake_log/{genes_to_work_on}_{sample}_samtools.log"
    benchmark: 
        snake_files + "/benchmark/{genes_to_work_on}_{sample}_samtools.bmk"
    shell:"""
        samtools view -bh -S {input.sam_transcriptome} > {output.bam_transcriptome} 2>> {log.snakemake_log}
        """

### Align to the whole transcriptome
rule transcriptome_alignment_whole:
    input:
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        sam_transcriptome = temp(base_path + "/results/alignments/transcriptome_whole/{sample}.sam"),
    params:
        transcriptome = transcriptome
    container:
        minimap2_img
    log:
        terminal_log = snake_files + "/log/{sample}_transcriptome_alignment.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_transcriptome_alignment.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_transcriptome_alignment.bmk"
    shell:"""
        minimap2 -ax map-ont --sam-hit-only \
        {params.transcriptome} \
        {input.fastq_filtered} > {output.sam_transcriptome} 2>> {log.snakemake_log}
        """

rule sam_to_bam_trans_whole:
    input:
        sam_transcriptome = base_path + "/results/alignments/transcriptome_whole/{sample}.sam"
    output:
        bam_transcriptome = base_path + "/results/alignments/transcriptome_whole/{sample}.bam"
    container:
        samtools_img
    log:
        terminal_log = snake_files + "/log/{sample}_samtools.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_samtools.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_samtools.bmk"
    shell:"""
        samtools view -bh -S {input.sam_transcriptome} > {output.bam_transcriptome} 2>> {log.snakemake_log}
        """

rule trigger_transcriptome_align:
    input:
        expand(base_path + "/results/alignments/transcriptome/{sample}_{genes_to_work_on}.bam",
        base_path = base_path, sample = samples, genes_to_work_on = genes_to_work_on),
        expand(base_path + "/results/alignments/transcriptome_whole/{sample}.bam",
        base_path = base_path, sample = samples)

