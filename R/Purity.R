#' Clustering Purity
#' 
#' Calculate the purity of the clustering results. For example, see 
#' \insertCite{Schaeffer_etal_2016_trust;textual}{funtimes}.
#' 
#' @details Following \insertCite{Manning_etal_2008;textual}{funtimes}, 
#' each cluster is assigned to the class which is most frequent in the cluster, then
#' \deqn{Purity(\Omega,C) = \frac{1}{N}\sum_{k}\max_{j}|\omega_k\cap c_j|,}
#' where  \eqn{\Omega=\{\omega_1,\ldots,\omega_K \}} is the set of identified 
#' clusters and \eqn{C=\{c_1,\ldots,c_J\}} is the set of classes. That is, within 
#' each class \eqn{j=1,\ldots,J} find the size of the most populous cluster from 
#' the \eqn{K-j} unassigned clusters. Then, sum together the \eqn{\min(K,J)} sizes 
#' found and divide by \eqn{N}, 
#' where \eqn{N} = \code{length(classes)} = \code{length(clusters)}.
#' 
#' If \eqn{\max_{j}|\omega_k\cap c_j|} is not unique for some \eqn{j}, 
#' it is assigned to the class which the second maximum is the smallest, to 
#' maximize the \eqn{Purity} (see `Examples').
#' 
#' The number of unique elements 
#' in \code{classes} and \code{clusters} may differ.
#' 
#' 
#' @param classes a vector with labels of true classes.
#' @param clusters a vector with labels of assigned clusters for which purity is to 
#' be tested. Should be of the same length as \code{classes}.
#' 
#' 
#' @return A list with two elements:
#' \item{pur}{purity value.}
#' \item{out}{table with \eqn{\min(K,J)} = \code{min(length(unique(classes)), 
#' length(unique(clusters)))} rows and the following columns: 
#' \code{ClassLabels}, \code{ClusterLabels}, and \code{ClusterSize}.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @keywords cluster
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' # Fix seed for reproducible simulations:
#' # RNGkind(sample.kind = "Rounding") #run this line to have same seed across R versions > R 3.6.0
#' set.seed(1)
#' 
#' ##### Example 1
#' #Create some classes and cluster labels:
#' classes <- rep(LETTERS[1:3], each = 5)
#' clusters <- sample(letters[1:5], length(classes), replace = TRUE)
#' 
#' #From the table below:
#' # - cluster 'b' corresponds to class A;
#' # - either of the clusters 'd' and 'e' can correspond to class B,
#' #   however, 'e' should be chosen, because cluster 'd' also highly 
#' #   intersects with Class C. Thus,
#' # - cluster 'd' corresponds to class C.
#' table(classes, clusters)
#' ##       clusters
#' ##classes a b c d e
#' ##      A 0 3 1 0 1
#' ##      B 1 0 0 2 2
#' ##      C 1 2 0 2 0
#' 
#' #The function does this choice automatically:
#' purity(classes, clusters)
#' 
#' #Sample output:
#' ##$pur
#' ##[1] 0.4666667
#' ##
#' ##$out
#' ##  ClassLabels ClusterLabels ClusterSize
#' ##1           A             b           3
#' ##2           B             e           2
#' ##3           C             d           2
#' 
#' 
#' ##### Example 2
#' #The labels can be also numeric:
#' classes <- rep(1:5, each = 3)
#' clusters <- sample(1:3, length(classes), replace = TRUE)
#' purity(classes, clusters)
#' 
.resolve_purity_tie <- function(contingency_table, max_indices) {
    # Helper function to resolve ties when multiple pairs have the same max frequency.
    # It chooses the pair where the cluster has the smallest second-most-frequent class,
    # thereby maximizing the chances for other classes to be matched well.
    if (any(dim(contingency_table) < 2)) {
        return(max_indices[1, , drop = FALSE])
    }
    
    # For each potential tied cluster, find the value of the second-most frequent class.
    second_max_values <- sapply(1:nrow(max_indices), function(j) {
        cluster_col_index <- max_indices[j, 2]
        sorted_col <- sort(contingency_table[, cluster_col_index], decreasing = TRUE)
        # If the column has fewer than 2 non-zero values, the second max is 0.
        if (length(sorted_col) < 2) return(0)
        return(sorted_col[2])
    })
    
    # Choose the tie that has the smallest second-max value.
    best_tie_index <- which.min(second_max_values)
    return(max_indices[best_tie_index, , drop = FALSE])
}

purity <- function(classes, clusters) 
{
  if (is.null(classes) || is.null(clusters))
    stop("classes and clusters must be provided.")
  if (!is.vector(classes) || !is.vector(clusters))
    stop("classes and clusters must be vectors.")
  if (length(classes) != length(clusters))
    stop("classes and clusters must have the same length.")
  if (length(classes) == 0L)
    stop("classes and clusters must be non-empty.")
  if (any(is.na(classes)) || any(is.na(clusters)))
    stop("classes and clusters must not contain missing values.")

  ClassLabels <- ClusterLabels <- ClusterSize <- numeric()
  contingency_table <- table(classes, clusters)
  i <- 1
  
  # Iteratively find the best class-cluster mapping
  while (min(dim(contingency_table)) > 0) {
    max_val <- max(contingency_table)
    max_indices <- which(contingency_table == max_val, arr.ind = TRUE)
    
    if (nrow(max_indices) > 1) { # Resolve ties
        best_index <- .resolve_purity_tie(contingency_table, max_indices)
    } else {
        best_index <- max_indices
    }
    
    # Assign the best pair
    ClassLabels[i] <- rownames(contingency_table)[best_index[1]]
    ClusterLabels[i] <- colnames(contingency_table)[best_index[2]]
    ClusterSize[i] <- max_val
    
    # Remove the assigned class and cluster from the table for the next iteration
    contingency_table <- contingency_table[
        -which(rownames(contingency_table) == ClassLabels[i]), 
        -which(colnames(contingency_table) == ClusterLabels[i]), 
        drop = FALSE
    ]
    
    # Re-evaluate the remaining table to remove all-zero rows or columns, if any
    contingency_table <- contingency_table[rowSums(contingency_table) > 0, , drop = FALSE]
    contingency_table <- contingency_table[, colSums(contingency_table) > 0, drop = FALSE]
    
    i <- i + 1
  }
  
  out <- data.frame(ClassLabels, ClusterLabels, ClusterSize)
  out <- out[order(out$ClassLabels),]
  list(pur = sum(ClusterSize)/length(classes), out = out)
}
