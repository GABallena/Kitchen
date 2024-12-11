import numpy as np
from scipy.stats import entropy

def calculate_conservation(alignment):
    # Conservation is the proportion of the most common residue in each column
    conservation_scores = [
        max(col.count(base) for base in set(col)) / len(col) for col in zip(*alignment)
    ]
    return np.mean(conservation_scores)

def calculate_entropy(alignment):
    # Calculate Shannon entropy for each column
    entropy_scores = [
        entropy([col.count(base) / len(col) for base in set(col)], base=2)
        for col in zip(*alignment)
    ]
    return np.mean(entropy_scores)

# Statistical model to combine metrics
def model_predict(conservation, entropy):
    return conservation - entropy  # Example: maximize conservation, minimize entropy

# Load alignment file
def load_alignment(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if not line.startswith(">")]

alignment = load_alignment("input.fasta")
conservation = calculate_conservation(alignment)
entropy_score = calculate_entropy(alignment)
quality_score = model_predict(conservation, entropy_score)

with open("output_statistical_quality_score.txt", "w") as f:
    f.write(f"Conservation: {conservation}\nEntropy: {entropy_score}\nQuality Score: {quality_score}\n")
