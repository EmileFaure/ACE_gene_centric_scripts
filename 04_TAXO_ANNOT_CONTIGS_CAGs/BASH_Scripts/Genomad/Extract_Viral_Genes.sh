cd /home/e/efaure/ACE/15_EXPLORE_GENOMAD

# Extract contigs ID (removing proviruses, viral sequences detected in host genomes)
cut -f1 /home/e/efaure/ACE/GENOMAD/genomad_output/all_fasta_concatenated_G_single_smaller_summary/all_fasta_concatenated_G_single_smaller_virus_summary.tsv | grep -v provirus > List_Contigs_Virus.tsv
tail -n+2 List_Contigs_Virus.tsv | sed 's/^/ACE_/' |  sed 's/$/_/' > tmp
mv tmp List_Contigs_Virus.tsv

# Extract same with taxo
cut -f1,11 /home/e/efaure/ACE/GENOMAD/genomad_output/all_fasta_concatenated_G_single_smaller_summary/all_fasta_concatenated_G_single_smaller_virus_summary.tsv | grep -v provirus | cut -f2 | tail -n+2 > List_Taxo_Virus.tsv

paste List_Contigs_Virus.tsv List_Taxo_Virus.tsv | sed 's/;/\t/g' > List_Contigs_Virus_withtaxo.tsv

# Add taxonomy of viruses to the full annotation table :
awk '
# Step 1: Read final_output.txt (store Contig ID → Taxonomy)
NR==FNR {
    tax[$1] = $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
    next
}

# Step 2: Process annotations.txt
FNR==1 {
    # Add taxonomy headers to the first row
    print $0 "\tGENOMAD_Domain\tGENOMAD_Realm\tGENOMAD_Kingdom\tGENOMAD_Phylum\tGENOMAD_Class\tGENOMAD_Order\tGENOMAD_Family"
    next
}

{
    # Extract contig ID from Gene ID (first two underscores define it)
    split($1, parts, "_")
    contig_id = parts[1] "_" parts[2] "_" parts[3] "_"  # Preserve the trailing underscore

    # Check if this contig ID exists in tax array (direct lookup, O(1))
    if (contig_id in tax) {
        print $0 "\t" tax[contig_id]  # Append taxonomy
    } else {
        print $0  # Print as is (no taxonomy)
    }
}' List_Contigs_Virus_withtaxo.tsv Annotation_Table_AGN_CDH_Tax_KEGG_EGG.tsv > Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD.tsv

# Extract ORFs / CD-Hit / AGC  with at least one viral gene
awk -F '\t' '$22!=""' Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD.tsv > Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_ViralORFsonly.tsv   # ORFs on viral contigs (1,848,248 ORFs)
cut -f2 Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_ViralORFsonly.tsv | tail -n+2 | sort | uniq > List_AGC_with_ViralORF.tsv #176,170 AGC
cut -f7 Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_ViralORFsonly.tsv | tail -n+2 | sort | uniq > List_CDHit_with_ViralORF.tsv #809,204 CD-Hit clusters

# We can now extract lines corresponding to viral AGC / CD-Hit from relative transformed matrices to sum them
head -n1 /home/e/efaure/ACE/10_APPLYTRANSFO/GM_CDHit_T60MAXGQ_relative-transformed.tsv > GM_VIRAL_CDHit_T60MAXGQ_relative-transformed.tsv
grep -Fwf List_CDHit_with_ViralORF.tsv /home/e/efaure/ACE/10_APPLYTRANSFO/GM_CDHit_T60MAXGQ_relative-transformed.tsv >> GM_VIRAL_CDHit_T60MAXGQ_relative-transformed.tsv

# Investigation of CAGs 29 and 79 :

sort list_contigs_CAG_29.tsv | uniq > list_unique_contigs_CAG_29.tsv
sort list_contigs_CAG_79.tsv | uniq > list_unique_contigs_CAG_79.tsv

grep -Fwf list_unique_contigs_CAG_29.tsv List_Contigs_Virus_withtaxo.tsv > list_contigs_Virus_withtaxo_CAG_29.tsv

grep -Fwf list_unique_contigs_CAG_79.tsv List_Contigs_Virus_withtaxo.tsv > list_contigs_Virus_withtaxo_CAG_79.tsv

# Focus on mimiviridae
grep Mimiviridae List_Contigs_Virus_withtaxo.tsv | cut -f1 > Contigs_Mimiviridae.tsv
grep -Fwf Contigs_Mimiviridae.tsv list_unique_contigs_CAG_79.tsv > Contigs_Mimiviridae_CAG_79.tsv
grep -Fwf Contigs_Mimiviridae.tsv list_unique_contigs_CAG_29.tsv > Contigs_Mimiviridae_CAG_29.tsv
cat Contigs_Mimiviridae_CAG_79.tsv Contigs_Mimiviridae_CAG_29.tsv | sort | uniq > Contigs_Mimiviridae_CAGs_Mertz.tsv

awk '
NR==FNR {ids[$1]=1; next}
{
    for (id in ids) {
        if (index($1, id) == 1) print
    }
}' Contigs_Mimiviridae_CAG_79.tsv Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_ViralORFsonly.tsv > Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_Mimiviridae_CAG79_only.tsv

# Focus on Autographiviridae

grep Autographiviridae List_Contigs_Virus_withtaxo.tsv | cut -f1 > Contigs_Autographiviridae.tsv
grep -Fwf Contigs_Autographiviridae.tsv list_unique_contigs_CAG_29.tsv > Contigs_Autographiviridae_CAG_29.tsv
grep -Fwf Contigs_Autographiviridae.tsv list_unique_contigs_CAG_79.tsv > Contigs_Autographiviridae_CAG_79.tsv
cat Contigs_Autographiviridae_CAG_29.tsv Contigs_Autographiviridae_CAG_79.tsv | sort | uniq > Contigs_Autographiviridae_CAG_Mertz.tsv

awk '
NR==FNR {ids[$1]=1; next}
{
    for (id in ids) {
        if (index($1, id) == 1) print
    }
}' Contigs_Autographiviridae_CAG_29.tsv Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_ViralORFsonly.tsv > Annotation_Table_AGN_CDH_Tax_KEGG_EGG_GENOMAD_Autographiviridae_CAG_29_only.tsv
