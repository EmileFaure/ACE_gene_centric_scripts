import pandas as pd
import csv
import ast

with open("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_ATT_RSquared15.txt","r") as data:
	dictionnary = ast.literal_eval(data.read())

with open("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_ATT_RSquared15.csv", 'w') as output:
	writer = csv.writer(output)
	for key, values in dictionnary.items():
		for value in values:
			writer.writerow([key,value])

