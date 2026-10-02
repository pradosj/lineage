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
library(lineage)
library(tidyverse)

# Build an gene expression matrix
counts <- t(as.matrix(iris[1:4]))
y <- iris$Species %>% fct_relevel("versicolor","virginica","setosa")

# Fit an linear ORM model in 2 pass
w <- orm_fit_pruned(counts,y,n=1,LAMBDA=1e-3)

pred <- orm_predict(counts,w)
tibble(pred,y) %>% ggplot(aes(x=y,y=pred)) + geom_boxplot()
enframe(w) %>% ggplot(aes(y=name,x=value)) + geom_col() + ggtitle("Linear model weights")
```


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


