perl -F'\t' -lane 'for (@F) { s/^\s+|\s+$//g }; print join("\t", @F)' contig_lineage.tsv > contig_lineage_clean.tsv

grep -Fwf list_contigs_CAG_79.tsv contig_lineage_clean.tsv > CAG_79_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_29.tsv contig_lineage_clean.tsv > CAG_29_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_83.tsv contig_lineage_clean.tsv > CAG_83_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_35.tsv contig_lineage_clean.tsv > CAG_35_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_22.tsv contig_lineage_clean.tsv > CAG_22_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_85270.tsv contig_lineage_clean.tsv > CAG_85270_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_137.tsv contig_lineage_clean.tsv > CAG_137_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_33.tsv contig_lineage_clean.tsv > CAG_33_GTDBKraken.tsv
grep -Fwf list_contigs_CAG_39.tsv contig_lineage_clean.tsv > CAG_39_GTDBKraken.tsv

cut -f6 CAG_79_GTDBKraken.tsv | sort | uniq -c | sort -n
#...

