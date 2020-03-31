# HAPLEX
Combinatorial and statistical prediction of gene expression from haplotype sequence

This repository contains the code for two methods that predict gene expression from haplotype data: HAPLEXR and HAPLEXD.
HAPLEXR is a penalized regression-based method that combines a SHAPEIT-like genetic model with allelic dosages to render continuous predictions.
HAPLEXD combines a suffix tree genetic model with spectral clustering to compute discrete expression predictions.


## HAPLEXR

To run the HAPLEXR example, clone the project, cd into the HAPLEXR directory, and run:

```bash
python HAPLEXR.py ../example_data/snp_annotations.txt ../example_data/gene_annotations.txt ../example_data/dosages.txt  ../example_data/clusterings.txt ../example_data/expression_phenotypes.txt example
```
## Logistic regression

To run a logistic regression example, clone the project, cd into the HAPLEX directory, and run:

```bash
Rscript LogisticRegression.R example_data/continuous/snp_annotations.txt example_data/continuous/gene_annotations.txt example_data/continuous/dosages.txt example_data/discrete/expression_discrete_E2.csv 2 0.5 example_data/example_predictions.txt
```
