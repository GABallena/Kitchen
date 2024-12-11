import numpy as np

def calculate_alignment_score(alignment):
    # Example scoring function: number of identical positions
    return sum(all(col[0] == c for c in col) for col in zip(*alignment))

def random_alignment_score(alignment):
    # Generate a random permutation of the alignment
    randomized = ["".join(np.random.permutation(seq)) for seq in alignment]
    return calculate_alignment_score(randomized)

def monte_carlo_alignment_test(alignment, n_simulations=1000):
    observed_score = calculate_alignment_score(alignment)
    simulated_scores = [random_alignment_score(alignment) for _ in range(n_simulations)]
    p_value = np.mean([sim >= observed_score for sim in simulated_scores])
    return observed_score, p_value

# Load alignment file
def load_alignment(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if not line.startswith(">")]

alignment = load_alignment("input.fasta")
score, p_value = monte_carlo_alignment_test(alignment)

with open("output_randomness_score.txt", "w") as f:
    f.write(f"Observed Score: {score}\nP-value: {p_value}\n")
