#' Clustering Purity
#' 
#' Calculate purity of the clustering results. For example, see 
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
#' it is assigned to the class which second maximum is the smallest, to 
#' maximize the \eqn{Purity} (see `Examples').
#' 
#' Number of unique elements 
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
purity <- function(classes, clusters) 
{
  ClassLabels <- ClusterLabels <- ClusterSize <- numeric()
  tmp <- table(classes, clusters)
  i <- 1
  while (min(dim(tmp)) > 0) {
    ind <- which(tmp == max(tmp), arr.ind = TRUE)
    if (nrow(ind) > 1) { #several equal max values
      if (any(dim(tmp) < 2)) {
        ind <- ind[1,]
      } else {
        #for each class with equal max, select the second max and choose minimum of that
        ind <- ind[which.min(sapply(1:nrow(ind), function(j) 
          sort(tmp[,ind[j,2]], decreasing = TRUE)[2])), ]
      }
    }
    ClassLabels[i] <- rownames(tmp)[ind[1]]
    ClusterLabels[i] <- colnames(tmp)[ind[2]]
    ClusterSize[i] <- max(tmp)
    tmp <- tmp[-which(rownames(tmp) == ClassLabels[i]), 
               -which(colnames(tmp) == ClusterLabels[i]), drop = FALSE]
    # Reevaluate the remaining table to remove all-zero rows or columns, if any
    tmp <- tmp[rowSums(tmp) != 0, , drop = FALSE]
    tmp <- tmp[, colSums(tmp) != 0, drop = FALSE]
    i <- i + 1
  }
  out <- data.frame(ClassLabels, ClusterLabels, ClusterSize)
  out <- out[order(out$ClassLabels),]
  list(pur = sum(ClusterSize)/length(classes), out = out)
}
