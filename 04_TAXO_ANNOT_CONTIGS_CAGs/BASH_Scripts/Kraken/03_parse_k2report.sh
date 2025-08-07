#!/bin/bash
awk -F'\t' '{print $7, $8}' all_fasta_concatenated_G_single_smaller.fa.k2report > K2report_lineage_map.txt
