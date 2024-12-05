import pandas as pd
from multiprocessing import Pool
from ann_linkage_clustering.lib import make_cags_with_ann
from ann_linkage_clustering.lib import iteratively_refine_cags
from ann_linkage_clustering.lib import make_nmslib_index

df=pd.read_table('/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.tsv',index_col=[0]).transpose()

# Maximum distance threshold (use any value)
max_dist=0.2

# Distance metric (only 'cosine' is supported)
distance_metric="cosine"

# Multiprocessing pool (pick any number of threads, in this case `1`)
threads = 4
pool = Pool(threads)

# Linkage type (only `average` is fully supported)
linkage_type = "average"

# Make the ANN index
index = make_nmslib_index(df)

# Make the CAGs in the first round
cags = make_cags_with_ann(
    index,
    max_dist,
    df,
    pool,
    threads=threads,
    distance_metric=distance_metric,
    linkage_type=linkage_type
)

# Iteratively refine the CAGs (this is the part that is hardedcoded to 
# use average linkage clustering, while the step above could technically
# use any of `complete`, `single`, `average`, etc.)
iteratively_refine_cags(
    cags,
    df.copy(),
    max_dist,
    distance_metric=distance_metric,
    linkage_type=linkage_type,
    threads=threads
)

# open file for writing
f = open("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_FL_RSquared10.txt","w")

# write file
f.write( str(cags) )

# close file
f.close()
