# Snakemake pipeline for phylogenetic analysis of metagenomic data
from snakemake.utils import min_version, validate
import os
import sys
from typing import List, Set

# Set minimum Snakemake version
min_version("6.0.0")

# Configuration and schema validation
configfile: "configs/config.yaml"

# Pipeline directories
DIRS = {
    'INPUT': 'data/metagenomic_reads',
    'RESULTS': 'results',
    'LOGS': 'logs',
    'BENCHMARKS': 'benchmarks'
}

# Ensure required directories exist
for d in DIRS.values():
    os.makedirs(d, exist_ok=True)

# Input validation
SAMPLES = glob_wildcards(os.path.join(DIRS['INPUT'], "{sample}.fastq.gz")).sample
if not SAMPLES:
    sys.exit("Error: No input samples found in data/metagenomic_reads/")

def validate_busco_output(scg_dir: str) -> bool:
    """Validate BUSCO output directory structure."""
    return os.path.exists(scg_dir) and any(f.endswith(".faa") for f in os.listdir(scg_dir))

# Define the main rule with logging
rule all:
    input:
        "results/speciation_gene_trees.txt",
        expand("results/sh_test/{gene}_SH_test_results.txt", gene=get_scgs),
        expand("results/speciation_trees/{gene}_tree.nwk", gene=get_scgs)
    log:
        os.path.join(DIRS['LOGS'], "pipeline_completion.log")
    benchmark:
        os.path.join(DIRS['BENCHMARKS'], "pipeline_completion.tsv")

# Checkpoint to aggregate SCGs from BUSCO results
checkpoint aggregate_scgs:
    input:
        expand("results/scgs/{sample}_scgs.fasta", sample=SAMPLES)
    output:
        gene_list = "results/scgs/scg_list.txt"
    log:
        os.path.join(DIRS['LOGS'], "aggregate_scgs.log")
    benchmark:
        os.path.join(DIRS['BENCHMARKS'], "aggregate_scgs.tsv")
    run:
        def get_scg_set(sample: str) -> Set[str]:
            """Extract SCG set from BUSCO output for a sample."""
            scg_dir = f"results/scgs/{sample}_busco/single_copy_busco_sequences"
            if not validate_busco_output(scg_dir):
                raise ValueError(f"Invalid BUSCO output directory for sample {sample}")
            return {f.replace(".faa", "") for f in os.listdir(scg_dir) 
                   if f.endswith(".faa")}

        try:
            # Get SCGs that are present in all samples
            scg_sets = [get_scg_set(sample) for sample in SAMPLES]
            common_scgs = set.intersection(*scg_sets)

            if not common_scgs:
                raise ValueError("No common SCGs found across samples")

            # Write sorted SCGs to output file
            with open(output.gene_list, 'w') as f:
                for scg in sorted(common_scgs):
                    f.write(f"{scg}\n")
                    
            with open(log[0], 'w') as f:
                f.write(f"Successfully identified {len(common_scgs)} common SCGs\n")
                
        except Exception as e:
            with open(log[0], 'w') as f:
                f.write(f"Error: {str(e)}\n")
            raise

def get_scgs(wildcards) -> List[str]:
    """Get list of SCGs from checkpoint output."""
    with open(checkpoints.aggregate_scgs.get().output.gene_list) as f:
        return [line.strip() for line in f]

onsuccess:
    with open(os.path.join(DIRS['LOGS'], "pipeline_success.log"), "w") as f:
        f.write("Pipeline completed successfully\n")

onerror:
    with open(os.path.join(DIRS['LOGS'], "pipeline_error.log"), "w") as f:
        f.write("Pipeline failed\n")

def get_scgs(wildcards):
    with open(checkpoints.aggregate_scgs.get().output[0]) as f:
        return [line.strip() for line in f]

# 0. Assemble metagenomic reads into contigs using MetaSPAdes
rule metaspades_assembly:
    input:
        reads="data/metagenomic_reads/{sample}.fastq.gz"
    output:
        contigs="results/contigs/{sample}_contigs.fasta"
    log:
        "logs/metaspades/{sample}.log"
    conda:
        "envs/spades_env.yml"
    shell:
        "metaspades.py -s {input.reads} -o results/contigs/{wildcards.sample}_spades -t {threads} -m 32000 &> {log}"

# 0.1 Predict genes from assembled contigs
rule gene_prediction:
    input:
        contigs="results/contigs/{sample}_contigs.fasta"
    output:
        genes="results/genes/{sample}_predicted_genes.fasta"
    log:
        "logs/prodigal/{sample}.log"
    conda:
        "envs/metawrap_env.yml"
    shell:
        "prodigal -i {input.contigs} -a {output.genes} -p meta &> {log}"

# 0.2 Run BUSCO to identify Single-Copy Genes (SCGs)
rule busco_analysis:
    input:
        genes="results/genes/{sample}_predicted_genes.fasta"
    output:
        scgs="results/scgs/{sample}_scgs.fasta"
    log:
        "logs/busco/{sample}.log"
    conda:
        "envs/busco_env.yml"
    params:
        lineage=config["busco"]["lineage_dataset"],
        mode="protein"
    shell:
        "busco -i {input.genes} -o {wildcards.sample}_busco -m {params.mode} -l {params.lineage} --out_path results/scgs -c 4 &> {log}"

# 1. Align SCG sequences using T-Coffee, CLUSTALW, and MAFFT
rule individual_alignments:
    input:
        "results/scgs/{gene}.fasta"
    output:
        t_coffee="results/individual_alignments/{gene}_tcoffee_aligned.fasta",
        clustalw="results/individual_alignments/{gene}_clustalw_aligned.fasta",
        mafft="results/individual_alignments/{gene}_mafft_aligned.fasta"
    conda:
        "envs/t-coffee_env.yml"
    shell:
        """
        t_coffee -seq {input} -output fasta -outfile {output.t_coffee};
        clustalw -INFILE={input} -OUTFILE={output.clustalw};
        mafft --auto {input} > {output.mafft}
        """

# 2. Combine alignments using M-Coffee
rule mcoffee_meta_alignment:
    input:
        t_coffee="results/individual_alignments/{gene}_tcoffee_aligned.fasta",
        clustalw="results/individual_alignments/{gene}_clustalw_aligned.fasta",
        mafft="results/individual_alignments/{gene}_mafft_aligned.fasta"
    output:
        "results/prealigned/{gene}_mcoffee_aligned.fasta"
    conda:
        "envs/t-coffee_env.yml"
    shell:
        "t_coffee -seq {input.t_coffee},{input.clustalw},{input.mafft} -mode mcoffee -output fasta -outfile {output}"

# 3. Translate nucleotide sequences to amino acids using EMBOSS transeq
rule translate_sequences:
    input:
        "results/scgs/{gene}.fasta"
    output:
        "results/translated/{gene}_translated.fasta"
    conda:
        "envs/emboss_env.yml"
    shell:
        "transeq -sequence {input} -outseq {output}"

# 4. Trim alignments using ClipKit
rule trim_alignment:
    input:
        "results/prealigned/{gene}_mcoffee_aligned.fasta"
    output:
        "results/trimmed/{gene}_trimmed_alignment.fasta"
    conda:
        "envs/clipkit_env.yml"
    shell:
        "clipkit {input} -o {output}"

# 5. Generate tree using FastTree
rule generate_fasttree:
    input:
        "results/trimmed/{gene}_trimmed_alignment.fasta"
    output:
        "results/trees/{gene}_fasttree.nwk"
    conda:
        "envs/fasttree_env.yml"
    shell:
        "fasttree -nt {input} > {output}"

# 6. Refine alignments using TPMA
rule align_sequences:
    input:
        gene_alignment="results/trimmed/{gene}_trimmed_alignment.fasta",
        guide_tree="results/trees/{gene}_fasttree.nwk"
    output:
        "results/aligned/{gene}_tpma_aligned.fasta"
    conda:
        "envs/phylo_env.yml"
    shell:
        "TPMA -i {input.gene_alignment} -t {input.guide_tree} -o {output}"

# 7. Evaluate alignments using Monte Carlo-based randomness detection
rule evaluate_randomness:
    input:
        "results/aligned/{gene}_tpma_aligned.fasta"
    output:
        "results/evaluation/{gene}_randomness_score.txt"
    script:
        "scripts/evaluate_randomness.py"

# 8. Evaluate alignments using statistical quality scoring
rule evaluate_statistical_quality:
    input:
        "results/aligned/{gene}_tpma_aligned.fasta"
    output:
        "results/evaluation/{gene}_statistical_quality_score.txt"
    script:
        "scripts/evaluate_statistical_quality.py"

# 9. Evaluate alignments using mutual information
rule evaluate_mutual_information:
    input:
        "results/translated/{gene}_translated.fasta"
    output:
        "results/evaluation/{gene}_mutual_information_score.txt"
    script:
        "scripts/evaluate_mutual_information.py"

# 10. Evaluate alignments using deep learning-based scoring
rule evaluate_deep_learning:
    input:
        "results/translated/{gene}_translated.fasta"
    output:
        "results/evaluation/{gene}_deep_learning_score.txt"
    script:
        "scripts/evaluate_deep_learning.py"

# 11. Model testing using ModelTest-NG
rule model_testing:
    input:
        "results/aligned/{gene}_tpma_aligned.fasta"
    output:
        "results/model/{gene}_model.txt"
    conda:
        "envs/modeltest-ng_env.yml"
    shell:
        "modeltest-ng -i {input} -o results/model/{wildcards.gene}"

# 12. Parse best-fitting model
rule parse_best_model:
    input:
        "results/model/{gene}_model.txt"
    output:
        "results/model/{gene}_best_model.txt"
    script:
        "scripts/parse_best_model.py"

# 13. Build phylogenetic tree for each gene using the best model
rule build_gene_tree:
    input:
        alignment="results/aligned/{gene}_tpma_aligned.fasta",
        model="results/model/{gene}_best_model.txt"
    output:
        "results/trees/{gene}_tree.nwk"
    conda:
        "envs/iqtree_env.yml"
    shell:
        "iqtree -s {input.alignment} -m $(cat {input.model}) -nt AUTO -pre results/trees/{wildcards.gene}_tree"

# 14. Combine trees into a DensiTree visualization
rule construct_densitree:
    input:
        expand("results/trees/{gene}_tree.nwk", gene=get_scgs)
    output:
        "results/densitree/densitree.png"
    script:
        "scripts/construct_densitree.py"

# 15. Infer speciation and duplication events
rule infer_events:
    input:
        "results/trees/{gene}_tree.nwk"
    output:
        "results/events/{gene}_inferred_events.txt"
    script:
        "scripts/inferred_events.py"

# 16. SH test to compare gene trees with the DensiTree
rule run_sh_test:
    input:
        alignment="results/aligned/{gene}_tpma_aligned.fasta",
        trees="results/densitree/trees_to_compare.txt",
        model="results/model/{gene}_best_model.txt"
    output:
        "results/sh_test/{gene}_SH_test_results.txt"
    conda:
        "envs/iqtree_env.yml"
    shell:
        "iqtree -s {input.alignment} -z {input.trees} -m $(cat {input.model}) --test SH -nt AUTO -pre results/sh_test/{wildcards.gene}_SH_test_results"

# 17. Identify gene trees occurring during speciation events
rule filter_speciation_trees:
    input:
        events="results/events/{gene}_inferred_events.txt",
        sh_test="results/sh_test/{gene}_SH_test_results.txt"
    output:
        "results/speciation_gene_trees.txt"
    script:
        "scripts/filter_speciation_trees.py"

# 18. Keep Newick files of trees associated with speciation events
rule keep_speciation_trees:
    input:
        speciation_list="results/speciation_gene_trees.txt",
        trees=expand("results/trees/{gene}_tree.nwk", gene=get_scgs)
    output:
        "results/speciation_trees/{gene}_tree.nwk"
    shell:
        "grep -Ff {input.speciation_list} {input.trees} > {output}"

# 19. Calculate diversity indices
rule calculate_diversity:
    input:
        expand(f"{bracken_output_dir}/{{sample}}.bracken", sample=SAMPLES)
    output:
        "diversity_matrices.tsv"
    conda:
        config["conda_env_diversity"]
    shell:
        """
        python scripts/calculate_diversity.py {bracken_output_dir} diversity_matrices.tsv
        """