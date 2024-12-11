import numpy as np
from sklearn.metrics import mutual_info_score

def calculate_mutual_information(alignment):
    # Calculate pairwise mutual information between columns
    columns = list(zip(*alignment))
    mi_scores = [
        mutual_info_score(columns[i], columns[j])
        for i in range(len(columns))
        for j in range(i + 1, len(columns))
    ]
    return np.mean(mi_scores)

# Load alignment file
def load_alignment(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if not line.startswith(">")]

alignment = load_alignment("input.fasta")
mi_score = calculate_mutual_information(alignment)

with open("output_mutual_information_score.txt", "w") as f:
    f.write(f"Mutual Information Score: {mi_score}\n")
