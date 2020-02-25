from tree import Tree
from util import *
from classify import *


def build_haplexD(trainingHap, trainingExp, testHap, k, n_clusters):
	t=Tree(trainingHap,trainingExp)
	nodes=t.leaves()
	predicted=classify(nodes, k, metric="chi2", testHap, n_clusters)
	return predicted