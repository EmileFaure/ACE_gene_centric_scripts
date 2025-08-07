#FeGenie results we'll investigate are given by contigs, with operons of Fe-related genes on each contigs given, separated by lines of ####

awk -F"_" '{print $2"_"$3"_"}' Contigs_SAR11GTDBUniRef_CAG22.ID > CAG22_togrepfegenie
awk -F"_" '{print $2"_"$3"_"}' Contigs_SAR11GTDBUniRef_CAG83.ID > CAG83_togrepfegenie
awk -F"_" '{print $2"_"$3"_"}' Contigs_SAR11GTDBUniRef_CAG35.ID > CAG35_togrepfegenie

grep -f CAG22_togrepfegenie FeGenie-geneSummary-clusters.csv > CAG22_contigs_FeGenie_clusters.csv

grep -f CAG83_togrepfegenie FeGenie-geneSummary-clusters.csv > CAG83_contigs_FeGenie_clusters.csv

grep -f CAG35_togrepfegenie FeGenie-geneSummary-clusters.csv > CAG35_contigs_FeGenie_clusters.csv

awk -F"," '{print $1}' CAG22_contigs_FeGenie_clusters.csv | sort | uniq -c | sort -nr | sed -E 's/^ *//; s/ /\t/' > CAG22_contigs_FeGenie_clusters.occurence

awk -F"," '{print $1}' CAG83_contigs_FeGenie_clusters.csv | sort | uniq -c | sort -nr | sed -E 's/^ *//; s/ /\t/' > CAG83_contigs_FeGenie_clusters.occurence

awk -F"," '{print $1}' CAG35_contigs_FeGenie_clusters.csv | sort | uniq -c | sort -nr | sed -E 's/^ *//; s/ /\t/' > CAG35_contigs_FeGenie_clusters.occurence

awk -F"," '{print $4}' CAG22_contigs_FeGenie_clusters.csv | sort | uniq -c | sort -nr | sed -E 's/^ *//; s/ /\t/' > CAG22_contigs_FeGenie_clusters.occurence.function

awk -F"," '{print $4}' CAG83_contigs_FeGenie_clusters.csv | sort | uniq -c | sort -nr | sed -E 's/^ *//; s/ /\t/' > CAG83_contigs_FeGenie_clusters.occurence.function

awk -F"," '{print $4}' CAG35_contigs_FeGenie_clusters.csv | sort | uniq -c | sort -nr | sed -E 's/^ *//; s/ /\t/' > CAG35_contigs_FeGenie_clusters.occurence.function

# number of contigs involved :
awk -F"," '{print $3}' CAG22_contigs_FeGenie_clusters.csv | awk -F"_" '{print $1"_"$2}' | sort | uniq | wc -l

awk -F"," '{print $3}' CAG35_contigs_FeGenie_clusters.csv | awk -F"_" '{print $1"_"$2}' | sort | uniq | wc -l

awk -F"," '{print $3}' CAG83_contigs_FeGenie_clusters.csv | awk -F"_" '{print $1"_"$2}' | sort | uniq | wc -l

