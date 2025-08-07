#!/usr/bin/env python3

import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Extract contig ID and taxid from Kraken2 output")
    parser.add_argument("-i", "--input", required=True, help="Kraken2 output file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV with contig_id and taxid")
    parser.add_argument("--use-names", action="store_true", help="Set if Kraken2 output used --use-names option")
    return parser.parse_args()

def extract_taxid(taxid_field, use_names):
    if use_names:
        match = re.search(r'taxid (\d+)', taxid_field)
        return match.group(1) if match else "NA"
    else:
        return taxid_field

def convert_kraken2(input_file, output_file, use_names):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write("contig_id\ttaxid\n")
        for line in infile:
            if line.startswith("C") or line.startswith("U"):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    contig_id = parts[1]
                    taxid = extract_taxid(parts[2], use_names)
                    outfile.write(f"{contig_id}\t{taxid}\n")

if __name__ == "__main__":
    args = parse_args()
    convert_kraken2(args.input, args.output, args.use_names)

