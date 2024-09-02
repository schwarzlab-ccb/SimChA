import pandas as pd
from copy import copy
import numpy as np
import random
from ete3 import Tree
import subprocess
from os.path import join


def generate_unique_node_names(tree):
    """
    Assigns unique names to internal nodes and leaves in the tree.
    """
    node_count = 0

    def traverse(node):
        nonlocal node_count
        if node.is_leaf():
            node.add_features(name=f"{node_count}")
            node_count += 1
        else:
            node.add_features(name=f"{node_count}")
            node_count += 1
            for child in node.children:
                traverse(child)

    traverse(tree)

def tree_to_edge_list(tree):
    edge_list = []

    def traverse(node):
        for child in node.children:
            edge_list.append((node.name, child.name, child.dist))
            traverse(child)

    traverse(tree)
    return edge_list

def edge_list_to_dataframe(edge_list):
    df = pd.DataFrame(edge_list, columns=['Parent', 'Leaf', 'Distance'])
    return df


def create_distance_matrix(tree):
    # Collect all nodes
    nodes = [node for node in tree.traverse()]
    node_names = [node.name for node in nodes]

    # Initialize distance matrix
    distance_matrix = pd.DataFrame(index=node_names, columns=node_names)

    # Fill the distance matrix
    for i, node1 in enumerate(nodes):
        for j, node2 in enumerate(nodes):
            if i <= j:  # Only compute for upper triangle (matrix is symmetric)
                distance = node1.get_distance(node2)
                distance_matrix.at[node1.name, node2.name] = distance
                distance_matrix.at[node2.name, node1.name] = distance

    return distance_matrix




def generate_random_binary_tree(num_leaves, distance_max, only_return_dataframe = True, print_tree = False):
    tree = Tree()
    tree.populate(num_leaves)
    
    for node in tree.traverse("postorder"): 
        if node.children:
            distance = random.randint(1, distance_max)
            for child in node.children:
                child.dist = distance
    
    edge_list_tree = copy(tree)
    generate_unique_node_names(edge_list_tree)
    df = edge_list_to_dataframe(tree_to_edge_list(edge_list_tree))
    
    if print_tree:
        print(tree)  
        print()  

    if only_return_dataframe:
        return df
    else:
        distance_matrix = create_distance_matrix(tree)

        return df, tree, distance_matrix

def generate_tree(output, i):
    # Create the subfolder for the outputs
    subdir = join(output, str(i+1))
    subprocess.run([f"mkdir -p {subdir}"], shell=True)
    # Generate the random binary tree for input into SimChA
    df = generate_random_binary_tree(num_leaves = 16, distance_max = 10, print_tree = False)
    # Rename the column headers
    df.rename(columns={'Parent': 'ParentID', 'Leaf': "ID"}, inplace=True)
    # SimChA needs the ID first, then ParentID
    df = df[["ID", "ParentID", "Distance"]]
    first_row = pd.DataFrame({"ID":[0], "ParentID": [-1], "Distance": [0]})
    df = pd.concat([first_row, df]).reset_index(drop=True)
    df["Distance"] = df["Distance"].astype(int)
    #print(df)
    df.to_csv(join(subdir, "tree.tsv"), sep="\t", index=False)
    return

def run_simcha(output, i):
    subdir = join(output, str(i+1))
    cmd = f"dotnet run --no-build --project SimChA -T {subdir}/tree.tsv -O {subdir} -C tom_inferred_noWGD.json"
    subprocess.run(cmd, shell=True)
    return

def tree_events():
    # Create the output folder
    output = "out_trees"
    subprocess.run([f"mkdir -p {output}"], shell=True)
    n=1000
    for i in range(n):
        generate_tree(output, i)
        run_simcha(output, i)
    subprocess.run(f"tar -zcvf {output}.tar.gz {output}", shell=True)
    return


if __name__ == "__main__":
    tree_events()

    


