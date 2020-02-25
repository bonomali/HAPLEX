import pandas as pd 
import numpy as np
from sklearn.cluster import spectral_clustering
from util import *


def classify(STSdf,k, metric, haplotypes, n_clusters):
	"""
	Takes in:
	STS: a pandas dataframe of signiticant tracts
	k: Number of "top" STS to look at
	metric: conditional entropy or chi2, 
	determines which metric to look at when building adjacency matrix
	haplotypes: haplotypes to predict the expression of. 

	returns: predictions
	"""
	affinity_matrix=np.zeros((len(haplotypes), len(haplotypes)))

	tractized_haplotypes=[]
	for h in haplotypes:
		tractized_haplotypes.append(tractize(h))

	STSdf.sort_values(by=metric, ascending=False, inplace=True, ignore_index=True)
	STSdf=STSdf.head(n=k)
	STS=STSdf["path"].tolist()

	for i in haplotypes:
		for j in haplotypes:
			count=0
			for s in STS:
				if (s in i) and (s in j):
					count+=1
			affinity_matrix[i][j]=count/len(STS)

	predicted_labels=spectral_clustering(affinity_matrix, n_clusters=n_clusters)




