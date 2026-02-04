##############
### FLAME ###
##############

### align 
rule minimap_for_flame:
    input:
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        bam_flame_sorted = base_path + "/results/flame/{sample}_sorted.bam",
        bai_flame_sorted = base_path + "/results/flame/{sample}_sorted.bam.bai"
    params:
        genome = genome
    threads:
        4
    log:
        terminal_log = snake_files + "/log/{sample}_flame_alignment.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flame_alignment.log"
    benchmark:
        snake_files + "/benchmark/{sample}_flame_alignment.bmk"
    shell:"""
    minimap2 -t {threads} -ax splice -k14 --secondary=no --splice-flank=no {params.genome} {input.fastq_filtered} | \
    samtools view -bS -@ 4 -m 2G - | \
    samtools sort -@ {threads} -o {output.bam_flame_sorted}
    samtools index {output.bam_flame_sorted} {output.bai_flame_sorted}
    """

### flame
rule flame:
    input:
        bam_flame_sorted = base_path + "/results/flame/{sample}_sorted.bam",
        bai_flame_sorted = base_path + "/results/flame/{sample}_sorted.bam.bai",
        gff3_filtered = base_path + "/ref_files/filtered_MANE.gff3"
    output:
        FLAME_done = base_path + "/results/flame/{sample}/FLAME_done.txt"
    params:
        flame_script = flame_script, ### OJO
        genome        = genome,
        gff3_filtered = base_path + "/ref_files/filtered.gff3",
        out_dir       = base_path + "/results/flame/{sample}",
        flame_config = flame_config, ### OJO
        minimap2_dir  = "/imppc/labs/eclab/ijarne/miniconda3/envs/FLAMES/bin/", ### OJO - NO CAL CANVIAR RES
        input_fasta = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz", 
        moved_input_fasta = base_path + "/results/flame/{sample}/merged.fastq.gz"
    threads: 
        10
    conda:
        flame_env
    log:
        terminal_log = snake_files + "/log/{sample}_flame.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flame.log"
    benchmark:
        snake_files + "/Run/benchmark/{sample}_flame.bmk"
    shell:"""
    ln -sf {params.input_fasta} {params.moved_input_fasta}
    
    cd {params.flame_script}

    ./python/bulk_long_pipeline.py \
    -a {input.gff3_filtered} \
    --genomefa {params.genome} \
    --outdir {params.out_dir} \
    --config_file {params.flame_config} \
    --minimap2_dir {params.minimap2_dir} \
    --inbam {input.bam_flame_sorted}

    echo "FLAME has been performed on sample {wildcards.sample}" > {output.FLAME_done}
    """

rule trigger_flame:
    input:
        expand(base_path + "/results/flame/{sample}/FLAME_done.txt",
        base_path = base_path, sample = samples),
        # expand(base_path + "/results/flame/{sample}_sorted.bam",
        # base_path = base_path, sample = samples),
        # expand(base_path + "/results/flame/{sample}_sorted.bam.bai",
        # base_path = base_path, sample = samples)

