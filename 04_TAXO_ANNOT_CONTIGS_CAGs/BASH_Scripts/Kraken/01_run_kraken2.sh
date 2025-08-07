sample="all_fasta_concatenated_G_single_smaller.fa"
module load kraken2
kraken2 --db /scratch/work/dataref/k2_gtdb_220 --thread 32 --report $sample.k2report --unclassified-out $sample_unclassified.fasta --classified-out $sample_classified.fasta --report-minimizer-data /scratch/work/ACE_PROJECT/02_REF_FASTA/all/all_fasta_concatenated_G_single_smaller.fa > $sample.GTDB.kraken2
module unload kraken2
