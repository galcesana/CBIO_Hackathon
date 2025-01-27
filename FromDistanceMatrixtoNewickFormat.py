#!/usr/bin/env python3

import sys
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo


def build_nj_tree_from_distance_matrix(distance_file, newick_outfile):
    """
    Reads a tab-delimited *full NxN symmetric* distance matrix file with format:
         (tab)  sample1  sample2  sample3 ...
    sample1     0        2        3       ...
    sample2     2        0        1       ...
    sample3     3        1        0       ...
    ...

    Converts the NxN matrix to the lower-triangle format Biopython expects,
    constructs a Neighbor-Joining tree, and writes the tree in Newick format
    to `newick_outfile`.
    """

    # Read lines, ignoring empty ones
    with open(distance_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # First line is the header row, e.g.:
    #   ['', 'sample1', 'sample2', 'sample3', ...]
    header_fields = lines[0].split('\t')
    # sample names (exclude the first column which might be empty or a label)
    sample_names = header_fields
    num_samples = len(sample_names)

    # Parse the NxN distance matrix rows
    # Each subsequent line i should be: sampleName_i, d(i,1), d(i,2), ...
    full_matrix = []
    for i in range(1, num_samples + 1):
        row_fields = lines[i].split('\t')
        # row_fields[0] should match sample_names[i-1]
        distances_str = row_fields[1:]  # distance values
        distances_float = [float(x) for x in distances_str]
        full_matrix.append(distances_float)

    # Convert the NxN symmetric matrix to Biopython's *lower-triangle* format.
    # For the i-th row (0-based), we keep columns [0..i-1].
    # So row i will have exactly i distance values (distance to previous samples).
    lower_triangle = []
    for i in range(num_samples):
        # Slices from the start of the row up to (but not including) i
        row_i = full_matrix[i][:i+1]
        lower_triangle.append(row_i)

    # Build a Biopython DistanceMatrix using the lower-triangle format
    dm = DistanceMatrix(names=sample_names, matrix=lower_triangle)

    # Construct a Neighbor-Joining tree
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)

    # Write tree in Newick format
    Phylo.write(nj_tree, newick_outfile, "newick")


def main(distance_file, newick_outfile):
    """
    Usage:
       python build_tree.py distance_matrix.tsv output_tree.nwk
    """
    # if len(sys.argv) != 3:
    #     print("Usage: python build_tree.py <distance_matrix.tsv> <output_tree.nwk>")
    #     sys.exit(1)

    # distance_file = sys.argv[1]
    # newick_outfile = sys.argv[2]

    build_nj_tree_from_distance_matrix(distance_file, newick_outfile)

    print(f"Neighbor-Joining tree written to: {newick_outfile}")


if __name__ == "__main__":
    main()