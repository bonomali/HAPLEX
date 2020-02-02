# HAPLEX
Combinatorial and statistical prediction of gene expression from haplotype sequence

This repository contains the code for two methods that predict gene expression from haplotype data: HAPLEXR and HAPLEXD.
HAPLEXR is a penalized regression-based method that combines a SHAPEIT-like genetic model with allelic dosages to render continuous predictions.
HAPLEXD combines a suffix tree genetic model with spectral clustering to compute discrete expression predictions.


## HAPLEXR

To run the HAPLEXR example, clone the project, cd into the HAPLEX directory, and run:

python HAPLEXR.py ../example_data/snp_annotations.txt ../example_data/gene_annotations.txt ../example_data/dosages.txt  ../example_data/clusterings.txt ../example_data/expression_phenotypes.txt example
