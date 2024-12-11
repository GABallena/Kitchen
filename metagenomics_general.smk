# Specify the configuration file
configfile: "configs/config.yaml"

import os
import re

# Ensure necessary directories exist
for directory in [
    "logs",
    config["trimmed_reads_dir"],
    config["fastqc_output_dir"],
    config["kraken_output_dir"],
    config["bracken_output_dir"],
    config["cleaning_results_dir"]
]:
    os.makedirs(directory, exist_ok=True)

# Load paths and configuration variables
raw_reads_dir = config["raw_reads_dir"]
trimmed_reads_dir = config["trimmed_reads_dir"]
fastqc_output_dir = config["fastqc_output_dir"]
kraken_output_dir = config["kraken_output_dir"]
bracken_output_dir = config["bracken_output_dir"]
kraken_db = config["kraken_db"]
cleaning_results_dir = config["cleaning_results_dir"]
adapter_file = config["adapter_file"]

# Dynamically collect sample names from raw_reads_dir
SAMPLES = []
for file in os.listdir(raw_reads_dir):
    if file.endswith('_1.fq.gz'):
        sample = re.match(r'(.+)_1\.fq\.gz', file)
        if sample:
            sample_name = sample.group(1)
            if os.path.exists(os.path.join(raw_reads_dir, f"{sample_name}_2.fq.gz")):
                SAMPLES.append(sample_name)


if not SAMPLES:
    raise ValueError(f"No valid raw read files found in '{raw_reads_dir}'.")

print(f"Found {len(SAMPLES)} valid sample pairs: {SAMPLES}")

rule all:
    input:
        # FastQC outputs
        expand(os.path.join("fastqc_reports", "{sample}_{read}_paired_fastqc.html"), sample=SAMPLES, read=["1", "2"]),
        
        # Trimmed reads
        expand(os.path.join(trimmed_reads_dir, "{sample}_{read}_paired.fastq.gz"), sample=SAMPLES, read=["1", "2"]),
        
        # Kraken2 outputs
        expand(os.path.join(kraken_output_dir, "{sample}.k2report"), sample=SAMPLES),
        expand(os.path.join(kraken_output_dir, "{sample}.kraken2"), sample=SAMPLES),
        
        # Bracken outputs
        expand(os.path.join(bracken_output_dir, "{sample}.bracken"), sample=SAMPLES),
        expand(os.path.join(bracken_output_dir, "{sample}.breport"), sample=SAMPLES)

rule run_trimming_cleaning_checking:
    input:
        R1=f"{config['raw_reads_dir']}/{{sample}}_1.fq.gz",
        R2=f"{config['raw_reads_dir']}/{{sample}}_2.fq.gz"
    output:
        paired_R1=temp(f"{config['trimmed_reads_dir']}/{{sample}}_1_paired.fastq"),
        paired_R2=temp(f"{config['trimmed_reads_dir']}/{{sample}}_2_paired.fastq"),
        unpaired_R1=temp(f"{config['trimmed_reads_dir']}/{{sample}}_1_unpaired.fastq"),
        unpaired_R2=temp(f"{config['trimmed_reads_dir']}/{{sample}}_2_unpaired.fastq"),
        final_R1=f"{config['trimmed_reads_dir']}/{{sample}}_1_paired.fastq.gz",
        final_R2=f"{config['trimmed_reads_dir']}/{{sample}}_2_paired.fastq.gz",
        qc_report=f"{config['cleaning_results_dir']}/{{sample}}_qc_report.txt"
    params:
        fastqc_dir=config["fastqc_output_dir"],
        adapter_file=config["adapter_file"],
        retries=config["retries"]  # Pass retries value
    conda:
        config["conda_env_trimming"]
    threads: config["threads_trimming"]  # Use threads from config.yaml
    log:
        "logs/trimming/{sample}.log"
    shell:
        """
        mkdir -p {config[trimmed_reads_dir]} {config[fastqc_output_dir]} {config[cleaning_results_dir]}
        python scripts/trimming_cleaning_checking.py \
            --input_R1 {input.R1} \
            --input_R2 {input.R2} \
            --output_paired1 {output.paired_R1} \
            --output_paired2 {output.paired_R2} \
            --unpaired1 {output.unpaired_R1} \
            --unpaired2 {output.unpaired_R2} \
            --fastqc_output_dir {params.fastqc_dir} \
            --adapter_file {params.adapter_file} \
            --sample {wildcards.sample} \
            --qc_report {output.qc_report} \
            --retries {params.retries} 2> {log}
        gzip < {output.paired_R1} > {output.final_R1}
        gzip < {output.paired_R2} > {output.final_R2}
        """


rule kraken_analysis:
    input:
        R1=f"{trimmed_reads_dir}/{{sample}}_1_paired.fastq.gz",
        R2=f"{trimmed_reads_dir}/{{sample}}_2_paired.fastq.gz"
    output:
        kraken_report=f"{kraken_output_dir}/{{sample}}.k2report",
        kraken_out=f"{kraken_output_dir}/{{sample}}.kraken2"
    params:
        db=config["kraken_db"],  # Pass the database path dynamically
    threads: config["threads_kraken"]  # Define threads directive dynamically
    resources:
        mem_mb=config["memory_kraken"]  # Memory in megabytes
    conda:
        config["conda_env_kraken"]
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --report {output.kraken_report} \
                --report-minimizer-data --minimum-hit-groups 3 \
                {input.R1} {input.R2} > {output.kraken_out}
        """

rule bracken_analysis:
    input:
        kraken_report=f"{kraken_output_dir}/{{sample}}.k2report"
    output:
        bracken_output=f"{bracken_output_dir}/{{sample}}.bracken",
        bracken_report=f"{bracken_output_dir}/{{sample}}.breport"
    conda:
        config["conda_env_kraken"]
    params: read_length=config["read_length_bracken"]
    shell:
        """
        bracken -d {kraken_db} -i {input.kraken_report} \
                -o {output.bracken_output} -w {output.bracken_report} \
                -r {params.read_length} -l S -t 1
        """
