#!/usr/bin/env python3

import sys

def vcf_to_treemut_matrix(vcf_file, alt_count_threshold=0):
    """
    Reads a multi-sample VCF (text) and returns a nested dict with
    presence/absence calls for each variant in each sample.

    :param vcf_file: Path to the VCF file (or '-' for stdin).
    :param alt_count_threshold: ALT count above which the variant is called 'present'.
    :return: (sample_list, variant_presence) where
             sample_list is a list of sample names in order
             variant_presence is dict of {variant_id: [0/1, 0/1, ...]} in the order of sample_list
    """

    # We will store sample names (in order) once we read the header line (#CHROM ...)
    sample_list = []

    # variant_presence will map variant_id -> list of 0/1 in the same sample order
    variant_presence = {}

    # Open file (or stdin) for reading
    if vcf_file == '-':
        f = sys.stdin
    else:
        f = open(vcf_file, 'r')

    for line in f:
        line=line.strip()

        # Skip meta-information lines
        if line.startswith('##'):
            continue

        # Parse the #CHROM header line to extract sample names
        if line.startswith('#CHROM'):
            # Example: #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1 sample2 ...
            header_fields = line.split('\t')
            # The first 9 columns are VCF fixed columns; the rest are samples
            sample_list = header_fields[9:]
            continue

        # Now we are in the body of the VCF
        fields = line.split('\t')
        chrom = fields[0]
        pos   = fields[1]
        var_id = fields[2]
        ref   = fields[3]
        alt   = fields[4]
        info  = fields[7]
        # Format is fields[8], then each sample column follows

        # We'll define a variant "identifier" that is unique per line
        # E.g. "chrom_pos_ref>alt" or something similar
        # You can choose your own labeling system for treemut input
        variant_identifier = f"{chrom}_{pos}_{ref}>{alt}"

        # The format column tells us the order of sub-fields, e.g. "REF_COUNT:ALT_COUNT:GENOTYPE"
        # We'll parse them accordingly.
        format_keys = fields[8].split(':')

        # We expect something like:
        #   format_keys = ["REF_COUNT", "ALT_COUNT", "GENOTYPE"]
        # for each sample. Let's figure out the index of REF_COUNT and ALT_COUNT
        try:
            ref_count_idx = format_keys.index("REF_COUNT")
            alt_count_idx = format_keys.index("ALT_COUNT")
        except ValueError as e:
            # If the format doesn't contain these fields, handle gracefully
            # or raise an error
            print("ERROR: REF_COUNT or ALT_COUNT not found in FORMAT field.")
            sys.exit(1)

        # Initialize a list of 0/1 calls for this variant
        presence_calls = []

        # For each sample column
        sample_calls = fields[9:]  # one element per sample

        for sample_call in sample_calls:
            # sample_call e.g. "16:0:0"
            subfields = sample_call.split(':')
            if len(subfields) <= alt_count_idx:
                # Malformed entry?
                presence_calls.append(0)
                continue

            try:
                alt_count = int(subfields[alt_count_idx])
            except ValueError:
                # If it's "?" or missing, we can default to absent or unknown
                alt_count = 0

            # presence if alt_count > alt_count_threshold
            if alt_count > alt_count_threshold:
                presence_calls.append(1)
            else:
                presence_calls.append(0)

        # Store in dictionary
        variant_presence[variant_identifier] = presence_calls

    if f is not sys.stdin:
        f.close()

    return sample_list, variant_presence


def main():
    """
    Command-line usage:
       python3 vcf_to_treemut_matrix.py input.vcf output_matrix.tsv
    """

    if len(sys.argv) < 3:
        print("Usage: python3 vcf_to_treemut_matrix.py <vcf_input> <matrix_output> [alt_count_threshold]")
        sys.exit(1)

    vcf_input = sys.argv[1]
    matrix_output = sys.argv[2]
    alt_count_threshold = 0
    if len(sys.argv) == 4:
        alt_count_threshold = int(sys.argv[3])

    sample_list, variant_presence = vcf_to_treemut_matrix(vcf_input, alt_count_threshold)

    # Write out the matrix as a tab-delimited file
    # First line:   variant_id  sample1   sample2   ...
    with open(matrix_output, 'w') as out:
        header = ["variant_id"] + sample_list
        out.write("\t".join(header) + "\n")

        for variant_id, presence_calls in variant_presence.items():
            row = [variant_id] + list(map(str, presence_calls))
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
