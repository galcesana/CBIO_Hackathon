# CBIO_Hackathon

Phylogenetic Tree Construction and Improvement
Overview
This project was developed as part of a hackathon for the CBIO course at the Hebrew University of Jerusalem. It provides a pipeline to generate and improve phylogenetic trees from VCF files by transforming them through a series of intermediate data representations.

Pipeline Description
The pipeline consists of the following main steps:

From VCF to Presence Matrix:

Converts a VCF file into a presence/absence matrix indicating the presence (1) or absence (0) of genetic variants for each sample.
Script: FromVCFtoPresenceMatrix.py
From Presence Matrix to Distance Matrix:

Computes pairwise distances (Hamming distances) between samples based on the presence/absence matrix.
Script: FromPresenceMatrixtoDistanceMatrix.py
From Distance Matrix to Newick Format:

Converts the distance matrix into a phylogenetic tree using the Neighbor-Joining (NJ) algorithm.
Script: FromDistanceMatrixtoNewickFormat.py
Improving the Phylogenetic Tree:

Enhances the tree by incorporating additional data from the presence matrix.
Script: newick_format_to_node_tree_1.py
Main Pipeline:

Automates the entire process for multiple samples.
Script: Main.py
File Descriptions
FromVCFtoPresenceMatrix.py:

Parses multi-sample VCF files to create a presence/absence matrix.
Allows thresholds for variant inclusion based on ALT counts.
FromPresenceMatrixtoDistanceMatrix.py:

Reads the presence matrix and generates a symmetric distance matrix using Hamming distance.
FromDistanceMatrixtoNewickFormat.py:

Reads a symmetric NxN distance matrix and constructs a Newick-formatted phylogenetic tree using the NJ method.
newick_format_to_node_tree_1.py:

Parses and processes Newick-formatted trees to integrate additional data for improvement.
Main.py:

Orchestrates the entire pipeline for a given list of samples.
Dependencies
The project requires the following Python libraries:

biopython
re
sys
Usage
To execute the pipeline for a list of samples:

Prepare a directory structure:

Place all VCF files in a vcf/ folder.
Create a Results/ folder for output files.
Run the main script:

bash
Copy
Edit
python Main.py
Output files will be generated in the Results/ folder:

Presence matrices (_\_presence.tsv)
Distance matrices (_\_distance.tsv)
Newick trees (_\_tree.nwk)
Improved Newick trees (_\_improved_tree.nwk)
Contributions
This project was developed collaboratively during a hackathon, showcasing teamwork and problem-solving in computational biology.
