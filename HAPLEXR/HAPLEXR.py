#!/usr/bin/env python
# coding: utf-8

# In[1]:


MODEL_NUM = "EXAMPLE"


# In[2]:


import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from scipy.stats import pearsonr
from multiprocessing import Pool, cpu_count
import sys
import os


# # Read input

# In[3]:
if len(sys.argv)==7:
    snp_file = sys.argv[1]
    gene_file = sys.argv[2]
    dosage_file = sys.argv[3]
    clusterings_file = sys.argv[4]
    expr_file = sys.argv[5]
    name=sys.argv[6]
else: 
    print("invalid input, need (in order) 6 arguments:")
    print("snp_annotations gene_annotations dosages clusterings expression name")
    print("see the example input files")
    exit(-1)

snp_annotations = pd.read_csv(snp_file)
gene_annotations = pd.read_csv(gene_file, index_col=0)
dosages = pd.read_csv(dosage_file, index_col=0)
clusterings = pd.read_csv(clusterings_file, index_col=0)
expression_phenotypes = pd.read_csv(expr_file, index_col=0)
tissue = name

# Differentiate haplotypes in index labels
clusterings.index = [start_var+"_h"+str(i%2) for i,start_var in enumerate(clusterings.index)]

# Split SNP annotations by chromosome for faster searching
snp_annotations_by_chr = {chr: snp_annotations[snp_annotations["chr"] == chr] for chr in set(snp_annotations["chr"])}


# # Model parameters

# In[4]:


SEED = 42

# Model parameters
ciswindow_size = 1e6
n_folds = 10
n_alphas = 10
tol = 1e-7


# # Define model and evaluation environment

# In[5]:


def get_haplotype_cluster_dummies(X):
    def _cat_to_vec(category, category_to_int_dict):
        categories = category_to_int_dict.keys()
        vec = [0] * (len(categories) - 1)
        category_num = category_to_int_dict[category]
        if category_num > 0:
            vec[category_num-1] = 1
        return np.array(vec)
    
    assert len(X[0]) == len(X[1]) and len(X)==2
    categories = list(set(list(X[0]) + list(X[1])))
    category_to_int_dict = {c:categories.index(c) for c in categories}
    
    dummies = []
    for i in range(len(X[0])):
        dummies.append(_cat_to_vec(X[0,i], category_to_int_dict) + _cat_to_vec(X[1,i], category_to_int_dict))
    return np.transpose(dummies)

def predict_holdout(gene):
    try:
        exppheno = expression_phenotypes.loc[gene].values

        # Select cis-snps and corresponding dosages
        gene_annotation = gene_annotations.loc[gene]
        gene_chr = gene_annotation["chr"]
        gene_start, gene_end = int(gene_annotation["start"]), int(gene_annotation["end"])
        ciswindow_start, ciswindow_end = gene_start - ciswindow_size, gene_end + ciswindow_size
        snp_annotations_chr = snp_annotations_by_chr[gene_chr]
        cissnp_annotations = snp_annotations_chr[snp_annotations_chr["pos"].between(ciswindow_start, ciswindow_end)]
        cissnps = list(cissnp_annotations["varID"])
        # Extract cissnp dosages
        cissnps_subset = list(set(cissnps).intersection(list(dosages.index)))
        cisgenos = dosages.loc[cissnps_subset,:] # TODO: Ensure all cissnps are in cisgenos (necessary?)
        # Extract cishaplo clusterings
        cissnps_haplo = []
        for cissnp in cissnps: cissnps_haplo += [cissnp+"_h0", cissnp+"_h1"]
        cissnps_haplo = [c for c in cissnps_haplo if c in clusterings.index]
        cisclusterings = clusterings.loc[cissnps_haplo,:]

        # Skip this gene if there are fewer than two variants
        if len(cisgenos) < 2:
            return [gene] + ["NA"]*(dosages.shape[1]//10), (gene, "NA", "NA", 0, 0, "CS<2")

        # Encode categorical clusterings with dummy variables
        for _ in range(0, len(cissnps_haplo), 2):    
            dummies = get_haplotype_cluster_dummies(cisclusterings.iloc[0:2,:].values)
            dummies_index = [cisclusterings.iloc[0:2,:].index[0].split("h0")[0]+"d"+str(i) for i in range(dummies.shape[0])]
            dummies_df = pd.DataFrame(dummies, index=dummies_index, columns=cisclusterings.columns)
            cisclusterings = cisclusterings.append(dummies_df)
            cisclusterings.drop(cisclusterings.index[0:2], axis=0, inplace=True)

        cisgenos = cisgenos.append(cisclusterings).T

        # Separate cisgenos and expphenos into train and test sets
        cisgenos_train = cisgenos.drop(cisgenos.index[::10], axis=0)
        cisgenos_test = cisgenos.loc[cisgenos.index[::10]]
        exppheno_train = [Y for i,Y in enumerate(exppheno) if i%10 != 0]
        exppheno_test = [Y for i,Y in enumerate(exppheno) if i%10 == 0]

        # Fit model for this gene
        model = ElasticNetCV(l1_ratio=0.5, n_alphas=n_alphas, cv=n_folds, tol=tol, random_state=SEED, n_jobs=1)
        model.fit(cisgenos_train.values, exppheno_train)

        if all(coef==0 for coef in model.coef_): # NOTE: Near-zero betas are not captured here
            r_train = "NA"
            return [gene] + ["NA"]*(dosages.shape[1]//10), (gene, "NA", "NA", 0, 0, "NB")
        
        # Save regression coefficients
        pd.DataFrame(model.coef_, index=cisgenos.columns).T.to_csv("output/model{}/betas/{}/{}.txt".format(MODEL_NUM, tissue, gene),
                                                                   sep=" ", index=False)
        n_dosage_nonzero_betas = np.count_nonzero(model.coef_[:-len(cisclusterings)])
        n_clustering_nonzero_betas = np.count_nonzero(model.coef_[-len(cisclusterings):])
        
        r_train = pearsonr(exppheno_train, model.predict(cisgenos_train.values))[0]

        exppheno_pred = model.predict(cisgenos_test.values)
        r_test = pearsonr(exppheno_test, exppheno_pred)[0]

        return [gene] + list(exppheno_pred), (gene, r_train, r_test, n_dosage_nonzero_betas, n_clustering_nonzero_betas, "OK")
    except:
        print("{}: unknown error".format(gene))
        return [gene] + ["NA"]*(dosages.shape[1]//10), (gene, "NA", "NA", "NA", "NA", "ERR")


# # Run model in parallel across genes

# In[6]:


os.makedirs("output/model{}/betas/{}/".format(MODEL_NUM, tissue), exist_ok=True)

genes = list(expression_phenotypes.index)
#with Pool(cpu_count()-2) as p:
#    model_outputs = p.map(predict_holdout, genes)

model_outputs = []
for gene in genes:
    model_outputs.append(predict_holdout(gene))


# # Save model predictions

# In[7]:


expression_phenotypes_pred_list = [o[0] for o in model_outputs]
expression_phenotypes_pred = pd.DataFrame(expression_phenotypes_pred_list, columns=["gene_id"] + list(dosages.columns[::10])).set_index("gene_id")
print(expression_phenotypes_pred.head())
fname = "output/model{}/{}_pred.txt".format(MODEL_NUM, tissue)
os.makedirs(os.path.dirname(fname), exist_ok=True)
expression_phenotypes_pred.to_csv(fname, sep='\t')
print('completed, see the output file: ', fname)

