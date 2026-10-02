# lineage

`lineage` is an R package to help infer lineage graph from single cell transcriptomic 
data where a pseudotime and a celltype score matrix is available. It also offers 
helper functions for pseudotime analysis with linear Ordinal Regression Model (ORM).



# Installation

``` r
devtools::install_github("pradosj/bmrm/bmrm")
devtools::install_github("BioinfoSupport/lineage/lineage")
```


# Usage

## Pseudotime analysis with linear Ordinal Regression Model (ORM)

``` r
w <- lineage::orm_fit_pruned(assay(trn,"counts"),trn$timepoint,n=25,LAMBDA=1e-3)
tst$orm_pred <- orm_predict(assay(tst,"counts"),w)
```

The above code shows how to fit an ORM model on a training set.
`LAMBDA` regularization parameters is a positive value controlling the "complexity" 
of the model and should be tuned for your data. Low `LAMBDA` tends to produce 
over-fitted models on the data, while high `LAMBDA` tends to produce 
non-specific model. 

After fitting, 25 genes with highest (resp. lowest) linear coefficients are kept 
(controlled by parameter `n`), and a second normalization and fitting pass 
is performed.



## Lineage with 

``` r
library(lineage)

# Generate a random pseudotime and identity score matrix
set.seed(123)
pseudotime <- runif(1000)
identity_scores <- matrix(runif(10000),1000,dimnames=list(NULL,LETTERS[1:10]))

# Find lineage
g <- lineage_graph_build(pseudotime,identity_scores) |>
  lineage_graph_prune()
plot_lineage_graph(g)
```


