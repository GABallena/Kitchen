import sys

def parse_events_file(events_file):
    """Parse the inferred events file."""
    speciation_nodes = set()
    with open(events_file) as f:
        for line in f:
            node, event = line.strip().split(": ")
            if event.lower() == "speciation":
                speciation_nodes.add(node)
    return speciation_nodes

def parse_sh_test_file(sh_test_file, significance_threshold=0.05):
    """Parse the SH test results file."""
    significant_trees = set()
    with open(sh_test_file) as f:
        for line in f:
            tree, p_value = line.strip().split(": ")
            if float(p_value) < significance_threshold:
                significant_trees.add(tree)
    return significant_trees

def filter_speciation_trees(events_file, sh_test_file, output_file):
    """Filter trees based on speciation events and SH test significance."""
    speciation_nodes = parse_events_file(events_file)
    significant_trees = parse_sh_test_file(sh_test_file)

    # Find the intersection
    valid_trees = speciation_nodes & significant_trees

    # Write the filtered trees to the output file
    with open(output_file, "w") as f:
        for tree in valid_trees:
            f.write(f"{tree}\n")

if __name__ == "__main__":
    # Command-line arguments
    events_file = sys.argv[1]  # e.g., "results/events/{gene}_inferred_events.txt"
    sh_test_file = sys.argv[2]  # e.g., "results/sh_test/{gene}_SH_test_results.txt"
    output_file = sys.argv[3]  # e.g., "results/speciation_gene_trees.txt"

    filter_speciation_trees(events_file, sh_test_file, output_file)
