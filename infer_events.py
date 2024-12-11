from ete3 import Tree

def infer_events(tree_file, output_file):
    """
    Infers orthology and paralogy relationships from a phylogenetic tree using the species-overlap algorithm.

    Args:
        tree_file (str): Path to the input phylogenetic tree file in Newick format.
        output_file (str): Path to save the inferred events.
    """
    # Load the phylogenetic tree
    tree = Tree(tree_file, format=1)

    results = []

    # Traverse the tree to infer events
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            # Get the set of species in the left and right partitions
            left_species = set(leaf.name.split("_")[0] for leaf in node.children[0].get_leaves())
            right_species = set(leaf.name.split("_")[0] for leaf in node.children[1].get_leaves())

            # Check for common species between the two partitions
            if left_species & right_species:
                event = "Duplication"
            else:
                event = "Speciation"

            # Record the event
            results.append({
                "node": node.name if node.name else "internal_node",
                "event": event,
                "left_species": left_species,
                "right_species": right_species,
            })

    # Save the results
    with open(output_file, "w") as f:
        f.write("Node\tEvent\tLeft_Species\tRight_Species\n")
        for res in results:
            f.write(f"{res['node']}\t{res['event']}\t{','.join(res['left_species'])}\t{','.join(res['right_species'])}\n")

    print(f"Inferred events saved to {output_file}.")

# Example usage
tree_file = "results/trees/sample_tree.nwk"
output_file = "results/inferred_events.txt"
infer_events(tree_file, output_file)
