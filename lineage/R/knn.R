
#' Select column elements in a matrix
#'
#' @param m matrix of values
#' @param i integer matrix with the same number of row as m, containing column-index
#'        elements to select in m
#'
#' @return matrix of with same dimension as i with values selected from m
#' @export
#' @examples
#' # Define a matrix of random value
#' m <- matrix(runif(1000),100,10)
#'
#' # Retrieve from each row the maximum value
#' row_select(m,max.col(m))
row_select <- function(m,i) {
	i <- as.matrix(i)
	array(m[cbind(as.vector(row(i)),as.vector(i))],dim=dim(i))
}


#' Order matrix rows
#'
#' @param m matrix whose rows have to be ordered.
#'
#' @return integer matrix with same dimension as m, with first column
#'  containing column-index in m with smallest value.
#' @export
#' @examples
#' # Define a matrix of random value
#' m <- matrix(runif(1000),100,10)
#'
#' # Order rows of the matrix
#' o <- row_order(m)
#'
#' # Retrieve sorted values
#' M <- row_select(m,o)
row_order <- function(m) {
	m <- as.matrix(m)
	o <- col(m)[order(as.vector(row(m)),as.vector(m))]
	matrix(o,byrow=TRUE,nrow=nrow(m),ncol=ncol(m))
}


#' Given a bin matrix, compute for each row number of element in each bin
#'
#' @param bin integer matrix of bins
#' @param nbins total number of bin
#'
#' @return matrix of size (n x nbins) with bin counts for each row
#' @export
row_tabulate <- function(bin,nbins=max(1,bin,na.rm=TRUE)) {
	bin <- as.matrix(bin)
	Matrix::sparseMatrix(
		i = row(bin),
		j = bin,
		x = rep_along(bin,1),
		dims = c(nrow(bin),nbins)
	)
}


#' Compute summary statistics for knn classifiers
#'
#' @param knn_index an integer matrix of size (n x k) containing indices of the
#'   k nearest neighbors
#' @param ref_labels a factor of labels for the reference dataset
#' @return a tibble with summary statistics
#' @export
#' @importFrom forcats as_factor
#' @importFrom tibble tibble
#' @examples
#' ref_coords <- matrix(runif(300*2),ncol=2)
#' ref_labels <- gl(3,100,labels = letters[1:3])
#' qry_coords <- matrix(runif(1000*2),ncol=2)
#' knn <- BiocNeighbors::queryKNN(ref_coords,qry_coords,k=25)
#' knn_summary_stats(knn$index,ref_labels)
knn_summary_stats <- function(knn_index,ref_labels) {
	ref_labels <- as_factor(ref_labels)
	n <- as.integer(ref_labels)[knn_index] |>
		array(dim(knn_index)) |>
		row_tabulate(nlevels(ref_labels)) |>
		as.matrix()
	colnames(n) <- levels(ref_labels)
	n <- n / ncol(knn_index)
	o <- row_order(-n)
	LAB <- factor(levels(ref_labels),levels(ref_labels))
	FRAC <- n[cbind(as.vector(row(n)),as.vector(o))] |> array(dim=dim(n))
	tibble(
		frac = n,
		order = o,
		label = LAB[o[,1]],
		label_frac = FRAC[,1],
		label2 = LAB[o[,2]],
		label2_frac = FRAC[,2]
	)
}

