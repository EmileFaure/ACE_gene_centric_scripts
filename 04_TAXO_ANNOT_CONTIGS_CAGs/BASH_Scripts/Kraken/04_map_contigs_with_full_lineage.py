#!/usr/bin/env python3

import argparse
import csv

def parse_lineage_map(lineage_file):
    taxid_to_lineage = {}
    lineage_stack = []

    with open(lineage_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue

            # Split the line at the first space only
            split_index = line.find(' ')
            if split_index == -1:
                continue  # skip malformed lines

            taxid = line[:split_index].strip()
            name_field = line[split_index + 1:]
            name = name_field.strip()

            # Compute depth based on indentation (2 spaces per level)
            depth = (len(name_field) - len(name)) // 2

            # Update lineage stack to correct depth
            lineage_stack = lineage_stack[:depth]
            lineage_stack.append(name)

            full_lineage = "; ".join(lineage_stack)
            taxid_to_lineage[taxid] = full_lineage

    return taxid_to_lineage

def map_contigs(contig_file, lineage_map_file, output_file):
    taxid_to_lineage = parse_lineage_map(lineage_map_file)

    with open(contig_file) as infile, open(output_file, 'w', newline='') as out:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["contig_id", "taxid", "lineage"])

        for row in reader:
            contig_id = row["contig_id"].strip()
            taxid = row["taxid"].strip()
            lineage = taxid_to_lineage.get(taxid, "NA")
            writer.writerow([contig_id, taxid, lineage])

    print(f"✅ Output written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Map contig taxids to full taxonomy lineages using a pre-parsed K2 report.")
    parser.add_argument("-c", "--contig-taxid", required=True, help="TSV file with contig_id and taxid")
    parser.add_argument("-l", "--lineage-map", required=True, help="Indented lineage map (K2report_lineage_map.txt)")
    parser.add_argument("-o", "--output", default="contig_lineage.tsv", help="Output file (TSV)")
    args = parser.parse_args()

    map_contigs(args.contig_taxid, args.lineage_map, args.output)

if __name__ == "__main__":
    main()

