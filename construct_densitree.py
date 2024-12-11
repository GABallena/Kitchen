from ete3 import Tree, TreeStyle

def construct_densitree(tree_files, output_image):
    trees = [Tree(tree_file) for tree_file in tree_files]

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.title.add_face(faces.TextFace("DensiTree Visualization", fsize=20), column=0)

    for t in trees:
        for node in t.traverse():
            nstyle = NodeStyle()
            nstyle["hz_line_color"] = "blue"
            nstyle["vt_line_color"] = "blue"
            nstyle["hz_line_width"] = 0.5
            nstyle["vt_line_width"] = 0.5
            node.set_style(nstyle)

    trees[0].render(output_image, tree_style=ts)

if __name__ == "__main__":
    import sys
    construct_densitree(sys.argv[1:], "results/densitree/densitree.png")
