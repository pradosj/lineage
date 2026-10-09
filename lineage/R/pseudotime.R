
#' @importFrom tibble tibble column_to_rownames
#' @importFrom dplyr mutate filter arrange left_join group_by cross_join count across if_else everything
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @importFrom tidyr pivot_wider
NULL

utils::globalVariables(c("group.x", "group.y", "y.x", "y.y", "n.x", "n.y"))

#' Compute a balanced cost matrix with optional grouping
#'
#' Grouping allow to train a single linear model performant on distinct set of samples like Male/Female
#'
#' @param y a factor of target labels
#' @param group a factor of sample groups whose length match length(y). It is used
#'        to update the cost matrix so that inter-group penalties vanish.
#' @export
#' @examples
#' orm_balanced_cost_matrix(iris$Species)
orm_balanced_cost_matrix <- function(y,group=NULL) {
	n <- tibble(y,group) %>%
		group_by(across(everything())) %>%
		dplyr::count() %>%
		mutate(lab = str_c(group,y,sep = "__")) %>%
		mutate(group = group %||% "all") %>%
		arrange(group,y)

	C <- cross_join(n,n) %>%
		filter(group.x==group.y) %>%
		mutate(n = if_else((group.x==group.y) & (as.integer(y.x)<as.integer(y.y)),n.x*n.y,0)) %>%
		mutate(loss = if_else(n>0,1/n,0)) %>%
		pivot_wider(id_cols="lab.x",names_from = "lab.y",values_from = "loss",values_fill=0) %>%
		column_to_rownames("lab.x") %>%
		as.matrix()
	C <- C*(sum(n$n)^2 - sum(n$n^2))/2

	list(
		C = C[n$lab,n$lab],
		y = factor(str_c(group,y,sep="__"),n$lab)
	)
}

#' Apply classical log-normalization on a raw counts
#'
#' @param counts the expression matrix, generally given as raw counts
#' @param scaling_factor the normlized L1-norm of the samples
#' @return log2(scale(counts,1/scaling_factor) + 1)
#' @export
#' @examples
#' counts <- t(as.matrix(iris[1:4]))
#' logcounts <- orm_lognorm_counts(counts)
#' colSums(2^logcounts - 1)
orm_lognorm_counts <- function(counts,scaling_factor=100) {
	log1p(scale(counts,center = FALSE,scale = pmax(colSums(counts),1)/scaling_factor))/log(2)
}


#' Normalize and return prediction of an ORM model for a given expression matrix
#'
#' The expression matrix is first reduced to genes in the model, and normlized
#'
#' @param counts gene x sample expression matrix with a maximum number of rownames
#'        matching with names(w)
#' @param w a named vector of the model linear weight
#' @param normalize_fun Normalization function to use before fitting
#' @return a numeric vector of predicted values for each column of counts
#' @export
#' @examples
#' orm_balanced_cost_matrix(iris$Species)
orm_predict <- function(counts,w,normalize_fun=orm_lognorm_counts) {
	stopifnot(is.numeric(w))
	stopifnot(!is.null(names(w)))
	w <- w[names(w) %in% rownames(counts)]
	as.vector(w %*% normalize_fun(counts[names(w),]))
}


#' Normalize an expression matrix and fit an ORM model
#'
#' @param counts gene x sample expression matrix
#' @param y A factor defining the ordinal labels for each sample in x
#' @param C Cost matrix to penalize miss-predictions between labels, C[i,j]
#'        being the cost for predicting label i instead of label j
#' @param normalize_fun Normalization function to use before fitting
#' @param ... additional parameters are passed to `orm_fit`
#' @return linear model weights
#' @export
#' @importFrom stats setNames
#' @examples
#' library(tidyverse)
#' counts <- t(as.matrix(iris[1:4]))
#' y <- iris$Species %>% fct_relevel("versicolor","virginica","setosa")
#' w <- orm_fit(counts,y,LAMBDA=1e-3)
#' pred <- orm_predict(counts,w)
#' tibble(pred,y) %>% ggplot(aes(x=y,y=pred)) + geom_boxplot()
orm_fit <- function(counts,y,C=NULL,normalize_fun=orm_lognorm_counts,...) {
	if (is.null(C)) {
		Cy <- orm_balanced_cost_matrix(y)
	} else {
		Cy <- list(C=C,y=y)
	}
	logcounts <- normalize_fun(counts)
	w <- bmrm::nrbm(bmrm::ordinalRegressionLoss(t(logcounts),y=Cy$y,C=Cy$C),...)
	w <- stats::setNames(as.vector(w),rownames(counts))
	w <- sort(w)
	w
}


#' Two pass fitting of a ORM model with pruning
#'
#' Fit a first ORM model with orm_fit() and select n top and n bottom features.
#' Then fit a second ORM model on selected features.
#' Note that in the second pass, the expression matrix is re-normlized only
#' considering selected features.
#'
#' @param counts gene x sample expression matrix
#' @param n Number of gene to keep after pruning (n positive + n negative)
#' @param ... additional parameters are passed to `orm_fit`, and should at least include y
#' @return linear model weights
#' @export
#' @importFrom utils head tail
#' @examples
#' library(tidyverse)
#' counts <- t(as.matrix(iris[1:4]))
#' y <- iris$Species %>% fct_relevel("versicolor","virginica","setosa")
#' w <- orm_fit_pruned(counts,y,n=1,LAMBDA=1e-3)
#' pred <- orm_predict(counts,w)
#' tibble(pred,y) %>% ggplot(aes(x=y,y=pred)) + geom_boxplot()
orm_fit_pruned <- function(counts,...,n=25L) {
	w <- orm_fit(counts,...)
	counts <- counts[c(names(utils::head(w,n=n)),names(utils::tail(w,n=n))),,drop=FALSE]
	w2 <- orm_fit(counts=counts,...)
	w2
}







