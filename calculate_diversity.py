import pandas as pd
import numpy as np
import yaml
import argparse
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import DistanceMatrix
from skbio.diversity.alpha import ace, margalef, menhinick, shannon, simpson, pielou_e
from skbio.diversity.beta import pw_distances, unweighted_unifrac, weighted_unifrac
from skbio.diversity import beta_diversity
from skbio.tree import TreeNode
from skbio import DistanceMatrix
from scipy.stats import entropy
from scipy.optimize import fsolve
import os
import glob

# Function to read .bracken files and extract species abundance data
def read_bracken(file):
    df = pd.read_csv(file, sep='\t')
    return df['new_est_reads'].values  # Only returning abundance

# Calculate alpha diversity using scikit-bio
def calculate_alpha_diversity(counts, metrics):
    results = {}
    for metric in metrics:
        if metric == 'ace':
            results[metric] = ace(counts)
        elif metric == 'margalef':
            results[metric] = margalef(counts)
        elif metric == 'menhinick':
            results[metric] = menhinick(counts)
        elif metric == 'shannon':
            results[metric] = shannon(counts)
        elif metric == 'simpson':
            results[metric] = simpson(counts)
        elif metric == 'pielou_e':
            results[metric] = pielou_e(counts)
        elif metric == 'inverse_simpson':
            results[metric] = 1 / simpson(counts)
        elif metric == 'gini_simpson':
            results[metric] = 1 - simpson(counts)
        elif metric == 'fisher_alpha':
            results[metric] = calculate_fisher_alpha(counts)
        elif metric == 'chao1':
            results[metric] = calculate_chao1(counts)
        elif metric == 'berger_parker':
            results[metric] = calculate_berger_parker(counts)
        else:
            alpha_values = alpha_diversity(metric, counts)
            results[metric] = alpha_values
    return results

# Custom implementation of Fisher's Alpha
def calculate_fisher_alpha(counts):
    counts = np.array(counts)
    def equation(a):
        return np.sum(np.log((a + np.arange(1, counts.sum()+1)) / a)) - counts.sum()
    
    initial_guess = 1.0
    fisher_alpha = fsolve(equation, initial_guess)[0]
    return fisher_alpha

# Custom implementation of Chao1 (using scipy or numpy)
def calculate_chao1(counts):
    counts = np.array(counts)
    singletons = np.sum(counts == 1)
    doubletons = np.sum(counts == 2)
    if doubletons > 0:
        chao1 = len(counts) + singletons**2 / (2 * doubletons)
    else:
        chao1 = len(counts) + singletons * (singletons - 1) / 2
    return chao1

# Custom implementation of Berger-Parker index
def calculate_berger_parker(counts):
    max_abundance = np.max(counts)
    total_abundance = np.sum(counts)
    return max_abundance / total_abundance

# Custom implementation of Sørensen-Dice Index
def calculate_sorensen_dice(counts1, counts2):
    intersection = np.sum(np.minimum(counts1, counts2))
    return 2 * intersection / (np.sum(counts1) + np.sum(counts2))

# Calculate beta diversity using scikit-bio
def calculate_beta_diversity(counts, metrics, tree=None):
    results = {}
    for metric in metrics:
        if metric in ['braycurtis', 'jaccard']:
            results[metric] = beta_diversity(metric, counts)
        elif metric == 'unweighted_unifrac' and tree:
            results[metric] = unweighted_unifrac(counts, tree)
        elif metric == 'weighted_unifrac' and tree:
            results[metric] = weighted_unifrac(counts, tree)
        elif metric == 'whittaker':
            results[metric] = whittaker_beta_diversity(counts)
    return results

# Custom implementation of Whittaker's Beta Diversity
def whittaker_beta_diversity(counts):
    total_species = np.sum(counts > 0, axis=1)
    gamma_diversity = np.sum(np.sum(counts, axis=0) > 0)
    alpha_diversity = np.mean(total_species)
    return (gamma_diversity - alpha_diversity) / alpha_diversity

# Custom implementation of Gamma Diversity
def calculate_gamma_diversity(counts):
    return np.sum(np.sum(counts, axis=0) > 0)

# Custom implementation of Faith's Phylogenetic Diversity (PD)
def calculate_faith_pd(counts, tree):
    return alpha_diversity('faith_pd', counts, tree=tree)

# Custom implementation of Mean Pairwise Distance (MPD)
def calculate_mpd(counts, tree):
    return alpha_diversity('mpd', counts, tree=tree)

# Custom implementation of Mean Nearest Taxon Distance (MNTD)
def calculate_mntd(counts, tree):
    return alpha_diversity('mntd', counts, tree=tree)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate diversity indices from bracken output and phylogenetic tree.")
    parser.add_argument("bracken_output_dir", help="Directory containing bracken output files")
    parser.add_argument("tree_files_dir", help="Directory containing phylogenetic tree files")
    parser.add_argument("output_file", help="Output file to save the diversity matrices")
    args = parser.parse_args()

    # Load configuration
    with open("config.yaml", "r") as config_file:
        config = yaml.safe_load(config_file)

    # Load sample files from the specified directory
    sample_files = glob.glob(os.path.join(args.bracken_output_dir, "*.bracken"))
    sample_names = [os.path.basename(file).replace(".bracken", "") for file in sample_files]

    # Load phylogenetic tree
    tree_files = glob.glob(os.path.join(args.tree_files_dir, "*.nwk"))
    if len(tree_files) != 1:
        raise ValueError("There should be exactly one Newick file in the tree_files_dir.")
    tree_file = tree_files[0]
    tree = TreeNode.read(tree_file) if os.path.exists(tree_file) else None

    # Read sample data
    sample_data = [read_bracken(file) for file in sample_files]

    # Pad sample data to ensure equal lengths
    max_length = max(len(sample) for sample in sample_data)
    padded_sample_data = [np.pad(sample, (0, max_length - len(sample)), 'constant') for sample in sample_data]

    # Alpha diversity metrics
    alpha_metrics = config["alpha_metrics"]
    alpha_results = {}

    # Calculate alpha diversity
    for i, data in enumerate(padded_sample_data):
        alpha_results[sample_names[i]] = calculate_alpha_diversity(data, alpha_metrics)

    # Prepare output for alpha diversity results
    alpha_output = "\nAlpha Diversity Results:\n"
    for sample, metrics in alpha_results.items():
        alpha_output += f"\nSample: {sample}\n"
        for metric, value in metrics.items():
            alpha_output += f"{metric}: {value}\n"

    # Beta diversity metrics
    beta_metrics = config["beta_metrics"]
    beta_results = calculate_beta_diversity(padded_sample_data, beta_metrics, tree)

    # Prepare output for beta diversity results
    beta_output = "\nBeta Diversity Results:\n"
    for metric, matrix in beta_results.items():
        beta_output += f"\n{metric.capitalize()} Distance Matrix:\n"
        distance_matrix = DistanceMatrix(matrix, ids=sample_names)
        beta_output += distance_matrix.to_data_frame().to_string() + "\n"

    # Calculate gamma diversity
    gamma_diversity = calculate_gamma_diversity(padded_sample_data)
    gamma_output = f"\nGamma Diversity: {gamma_diversity}\n"

    # Save results to a file
    with open(args.output_file, "w") as f:
        f.write(alpha_output + beta_output + gamma_output)

    # Print all results (alpha, beta, and gamma diversity) to stdout
    print(alpha_output + beta_output + gamma_output)

if __name__ == "__main__":
    main()