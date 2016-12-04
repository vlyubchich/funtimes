\name{Purity}
\alias{Purity}

\title{Clustering Purity}

\description{
Calculate purity of the clustering results.
}

\usage{
Purity(classes, clusters)
}

\arguments{
  \item{classes, clusters}{vectors of equal lengths, with labels of underlying true \code{classes} and assigned \code{clusters}. Number of unique elements in \code{classes} and \code{clusters} may differ.}
}

\details{
Following Manning et al. (2008), each cluster is assigned to the class which is most frequent in the cluster, then
\deqn{Purity(\Omega,C) = \frac{1}{N}\sum_{k}\max_{j}|\omega_k\cap c_j|,}
where  \eqn{\Omega=\{\omega_1,\ldots,\omega_K \}} is the set of identified clusters and \eqn{C=\{c_1,\ldots,c_J\}} is the set of classes. That is, within each class \eqn{j=1,\ldots,J} find the size of the most populous cluster from the \eqn{K-j} unassigned clusters. Then, sum together the \eqn{\min(K,J)} sizes found and divide by \eqn{N}, where \eqn{N} = \code{length(classes)} = \code{length(clusters)}.

If \eqn{\max_{j}|\omega_k\cap c_j|} is not unique for some \eqn{j}, it is assigned to the class which second maximum is the smallest, to maximize the \eqn{Purity} (see `Examples').
}

\value{
A list with two elements:
\item{purity}{purity value.}

\item{out}{table with \eqn{\min(K,J)} = \code{min(length(unique(classes)), length(unique(clusters)))} rows and the following columns: \code{ClassLabels}, \code{ClusterLabels}, and \code{ClusterSize}.}
}

\references{
Manning, C. D., Raghavan, P., and  Schutze, H. (2008). \emph{Introduction to
Information Retrieval}. New York: Cambridge University Press.

Schaeffer, E. D., Testa, J. M., Gel, Y. R., and Lyubchich, V. (2016). On information criteria for dynamic spatio-temporal clustering. In A. Banerjee et al. (Eds.) \emph{Proceedings of the 6th International Workshop on Climate Informatics: CI 2016}. NCAR Technical Note NCAR/TN-529+PROC, September 2016, pages 5--8. DOI: \href{http://dx.doi.org/10.5065/D6K072N6}{10.5065/D6K072N6}
}

\author{
Vyacheslav Lyubchich
}

\examples{
# Fix seed for reproduceable simulations:
set.seed(1)

##### Example 1
#Create some classes and cluster labels:
classes <- rep(LETTERS[1:3], each=5)
clusters <- sample(letters[1:5], length(classes), replace=TRUE)

#From the table below, it is clear that 
#for j=1 (Class A), cluster 'b' should be chosen, however,
#for j=2 (Class B), between clusters 'd' and 'e', 'e' should be chosen,
#because cluster 'd' also highly intersects with Class C:
table(classes, clusters)
##       clusters
##classes a b c d e
##      A 0 3 1 0 1
##      B 1 0 0 2 2
##      C 1 2 0 2 0

#The function does this choice automatically:
Purity(classes, clusters)

#Sample output:
##$purity
##[1] 0.4666667
##
##$out
##  ClassLabels ClusterLabels ClusterSize
##1           A             b           3
##2           B             e           2
##3           C             d           2


##### Example 2
#The labels can be also numeric:
classes <- rep(1:5, each=3)
clusters <- sample(1:3, length(classes), replace=TRUE)
Purity(classes, clusters)
}

\keyword{cluster}

