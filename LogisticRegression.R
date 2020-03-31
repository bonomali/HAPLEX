#!/usr/bin/env Rscript

library(glmnet)
set.seed(42)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  print("invalid input, need (in order) 7 arguments:")
  print("snp_annotations gene_annotations dosages expression name E alpha predictions_fn")
  quit(status=-1)
}
snp_file <- args[1]
gene_file <- args[2]
dosage_file <- args[3]
expr_file <- args[4]
n_expr_bins <- as.integer(args[5])
alpha <- as.numeric(args[6])
predictions_fn <- args[7]

expression <- as.matrix(read.csv(expr_file, row.names=1))
class(expression) <- "numeric"

# Select individuals with genotype and expression samples
genotype <- read.csv(dosage_file, row.names="Id", stringsAsFactors=F)
individuals <- intersect(colnames(expression), colnames(genotype))
genotype <- t(subset(genotype, select=individuals))
expression <- t(subset(expression, select=individuals))

# Select genes that have annotations
donors <- rownames(expression)
genes <- colnames(expression)
gene_annotations <- read.csv(gene_file, row.names=1)
expression <- expression[, intersect(genes, rownames(gene_annotations))]

snp_annotations <- read.csv(snp_file)

for (i in 1:length(genes)) {
  gene <- genes[i]
  
  # Select cis-SNP dosages
  geneinfo <- gene_annotations[gene,]
  start <- geneinfo$start - 1e6
  end <- geneinfo$end + 1e6
  cissnps <- subset(snp_annotations, snp_annotations$pos >= start & snp_annotations$pos <= end)
  cisgenos <- genotype[,intersect(colnames(genotype), cissnps$varID), drop = FALSE]
  
  # Separate individuals into train and test sets
  test_col_nums <- seq(1, nrow(cisgenos), 10)
  train_col_nums <- setdiff(1:nrow(cisgenos), test_col_nums)
  
  # Write the prediction file header if this is the first gene
  if (i==1)
    write(c("gene_id", individuals[test_col_nums]), file=predictions_fn, ncol=length(test_col_nums)+1, sep="\t")
  
  if (ncol(cisgenos) < 2) {
    print(paste(gene, "<2 cis-snps"))
    pred <- rep("NA", length(test_col_nums))
    write(c(gene, pred), predictions_fn, ncolumns=length(test_col_nums)+1, append=T, sep="\t")
    next
  }
  
  # Split dosages and expression into train and test sets
  exppheno <- expression[,gene]
  train_cisgenos <- cisgenos[train_col_nums,]
  train_exppheno <- exppheno[train_col_nums]
  test_cisgenos <- cisgenos[test_col_nums,]
  test_exppheno <- exppheno[test_col_nums]
  
  groupid <- sample(1:10, length(donors)-length(test_col_nums), replace=T)
  
  tryCatch({
    # Fit cross-validated elastic net to select lambda
    fit <- cv.glmnet(as.matrix(train_cisgenos), train_exppheno,
                      family = "multinomial", alpha = alpha,
                      nfolds = 10, foldid = groupid)
    nrow.best <- which.min(fit$cvm)
    lambda.best <- fit$lambda.min
    
    nonzero_betas <- FALSE
    for (k in 1:n_expr_bins) {
        if (any(fit$glmnet.fit$beta[[k]][,nrow.best] != 0.0))
        	 nonzero_betas <- TRUE
    }
    },
    error = function(cond) {
    	print(cond)
      pred <- rep("NA", length(test_col_nums))
  	  write(c(gene, pred), predictions_fn, ncolumns=length(test_col_nums)+1, append=T, sep="\t")
      next
    }
  )
  
  if (!nonzero_betas) {
    print(paste(gene, "all-0 betas"))
    pred <- rep("NA", length(test_col_nums))
    write(c(gene, pred), predictions_fn, ncolumns=length(test_col_nums)+1, append=T, sep="\t")
    next
  }
  
  # Predict on test set using the best lambda
  pred <- predict(fit, test_cisgenos, s=lambda.best, type="class")
  pred <- as.character(as.numeric(pred))
  write(c(gene, pred), predictions_fn, ncolumns=length(test_col_nums)+1, append=T, sep="\t")
}