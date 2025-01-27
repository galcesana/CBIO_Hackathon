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

    sample_list = []
    variant_presence = {}

    # Open file (or stdin) for reading
    if vcf_file == '-':
        f = sys.stdin
    else:
        f = open(vcf_file, 'r')

    for line in f:
        line = line.strip()

        # Skip meta-information lines
        if line.startswith('##'):
            continue

        # Parse the #CHROM header line to extract sample names
        if line.startswith('#CHROM'):
            header_fields = line.split('\t')
            # The first 9 columns are VCF fixed columns; the rest are samples
            sample_list = header_fields[9:]
            continue

        # Now we are in the body of the VCF
        fields = line.split('\t')
        chrom = fields[0]
        pos   = fields[1][3:]
        var_id = fields[2]
        ref   = fields[3]
        alt   = fields[4]
        info  = fields[7]


        if 'GENE=JAK2' in info:
            index = var_id[7:]
        # -------------------------------------------------------------
        # Skip the variant if IS_IN_EXCLUDED_REGION=1 is present in INFO
        if 'IS_IN_EXCLUDED_REGION=1' in info:
            continue
        # -------------------------------------------------------------
        # if 'MUT_TYPE=INDEL' in info:
        #     continue
        # if '6' in chrom:
        #     continue
        # if '7' in chrom:
        #     continue
        # if 'GENE=-' in info:
        #     continue
        # Create a variant identifier (you can adjust as you like)
        variant_identifier = f"{chrom}_{pos}_{ref}>{alt}"

        # The FORMAT column
        format_keys = fields[8].split(':')

        # Attempt to find REF_COUNT and ALT_COUNT in the format fields
        try:
            ref_count_idx = format_keys.index("REF_COUNT")
            alt_count_idx = format_keys.index("ALT_COUNT")
        except ValueError:
            print("ERROR: REF_COUNT or ALT_COUNT not found in FORMAT field.")
            sys.exit(1)

        presence_calls = []
        sample_calls = fields[9:]  # each element corresponds to a sample’s call

        for sample_call in sample_calls:
            subfields = sample_call.split(':')

            # If the subfields are malformed or missing, skip gracefully
            if len(subfields) <= alt_count_idx:
                presence_calls.append(0)
                continue

            try:
                alt_count = int(subfields[alt_count_idx])
            except ValueError:
                # If it's something like '?' or missing
                alt_count = 0

            if alt_count > alt_count_threshold:
                presence_calls.append(1)
            else:
                presence_calls.append(0)

        variant_presence[variant_identifier] = presence_calls

    if f is not sys.stdin:
        f.close()

    return sample_list, variant_presence


def main(vcf_input, matrix_output):
    """
    Command-line usage:
       python3 vcf_to_treemut_matrix.py input.vcf output_matrix.tsv
    """
    alt_count_threshold = 0
    # if len(sys.argv) == 4:
    #     alt_count_threshold = int(sys.argv[3])

    sample_list, variant_presence = vcf_to_treemut_matrix(vcf_input, alt_count_threshold)

    # Write out the matrix
    with open(matrix_output, 'w') as out:
        header = ["variant_id"] + sample_list
        out.write("\t".join(header) + "\n")

        for variant_id, presence_calls in variant_presence.items():
            row = [variant_id] + list(map(str, presence_calls))
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 vcf_to_treemut_matrix.py <vcf_input> <matrix_output> [alt_count_threshold]")
        sys.exit(1)

    vcf_input = sys.argv[1]
    matrix_output = sys.argv[2]

    main(vcf_input, matrix_output)
