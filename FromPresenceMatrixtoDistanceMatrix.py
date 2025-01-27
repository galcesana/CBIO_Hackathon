#!/usr/bin/env python3

import sys


def presence_matrix_to_distance_matrix(presence_matrix_file, output_file):
    """
    Reads a presence/absence matrix (tab-delimited) and computes a pairwise
    distance matrix among samples. The presence/absence matrix should have:

       1) The first row as headers:  variant_id  sample1  sample2  ...
       2) Subsequent rows: variant_id, followed by 0/1 calls for each sample.

    This function calculates Hamming distances, i.e., for each pair of samples,
    distance = sum of absolute differences in 0/1 calls over all variants.

    :param presence_matrix_file: Path to the tab-delimited presence matrix file
    :param output_file: Path to the output file where the distance matrix is written
    """

    # Read the presence/absence file
    with open(presence_matrix_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # First line is header: variant_id, sample1, sample2, ...
    header = lines[0].split('\t')
    sample_names = header[1:]  # Exclude 'variant_id'
    num_samples = len(sample_names)

    # Initialize a list (or dict) to hold each sample's presence/absence vector
    # Key: sample_name, Value: list of 0/1 for each variant
    sample_vectors = {sample: [] for sample in sample_names}

    # Read each subsequent line (one per variant)
    for line in lines[1:]:
        fields = line.split('\t')
        variant_id = fields[0]  # not needed for distance calculations, but present
        calls = fields[1:]  # the 0/1 calls for each sample

        # Store the presence/absence in each sample's vector
        for i, sample in enumerate(sample_names):
            sample_vectors[sample].append(int(calls[i]))

    # Now compute pairwise distances (Hamming distance) among samples
    # We'll create a matrix where distance_matrix[i][j] = distance(sample_i, sample_j)
    distance_matrix = [[0] * num_samples for _ in range(num_samples)]

    # Pre-convert each sample's presence calls to a list for fast iteration
    sample_call_lists = [sample_vectors[sample] for sample in sample_names]

    for i in range(num_samples):
        for j in range(i, num_samples):
            if i == j:
                distance = 0
            else:
                # Compare presence/absence vectors
                vec_i = sample_call_lists[i]
                vec_j = sample_call_lists[j]
                distance = sum(abs(a - b) for a, b in zip(vec_i, vec_j))

            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance  # symmetric

    # Write out the distance matrix
    # We'll write a tab-delimited matrix, with samples in rows and columns
    # First row = ['sample'] + sample_names
    with open(output_file, 'w') as out:
        # Write header
        out.write("\t" + "\t".join(sample_names) + "\n")
        # Write each row
        for i, sample_i in enumerate(sample_names):
            row_values = [str(distance_matrix[i][j]) for j in range(num_samples)]
            out.write(sample_i + "\t" + "\t".join(row_values) + "\n")


def main():
    """
    Command-line usage:
        python3 presence_matrix_to_distance_matrix.py presence_matrix.tsv distance_matrix.tsv
    """
    if len(sys.argv) != 3:
        print("Usage: python3 presence_matrix_to_distance_matrix.py <presence_matrix.tsv> <distance_matrix.tsv>")
        sys.exit(1)

    presence_matrix_file = sys.argv[1]
    output_file = sys.argv[2]

    presence_matrix_to_distance_matrix(presence_matrix_file, output_file)


if __name__ == "__main__":
    main()