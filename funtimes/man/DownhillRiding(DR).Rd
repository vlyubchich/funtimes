\name{DR}
\alias{DR}

\title{Downhill Riding (DR) Procedure}

\description{
Downhill Riding procedure for selecting optimal tuning parameters in clustering algorithms, using an (in)stability probe.
}
            
\usage{
DR(X, method, minPts = 3, theta = 0.9, B = 500, lb = -30, ub = 10)
}

\arguments{
\item{X}{a \eqn{n\times k} matrix where columns are \eqn{k} objects to be clustered, and each object contains n observations (objects could be a set of time series).}

\item{method}{the clustering method to be used (currently either \dQuote{TRUST} or \dQuote{DBSCAN}). If method is \code{DBSCAN}, then set \code{MinPts} and optimal \eqn{\epsilon} is selected using DR. If method is \code{TRUST}, then set \code{theta} and optimal \eqn{\delta} is selected using DR.}

\item{minPts}{the minimum number of samples in an \eqn{\epsilon}-neighborhood of a point to be considered as a core point. The \code{minPts} is to be used only with \code{DBSCAN} method. Default value is 3.}

\item{theta}{connectivity parameter \eqn{\theta \in (0,1)}, which is to be used only with \code{TRUST} method. Default value is 0.9.}

\item{B}{number of random splits in calculating the \code{Average Cluster Deviation (ACD)}. Default value is 500.}

\item{lb, ub}{end points for a range of search for the optimal parameter.}
}
              
\value{
A list containing the following components:

\item{P_opt}{the value of optimal parameter. If method is \code{DBSCAN}, then \code{P_opt} is optimal \eqn{\epsilon}. 
If method is \code{TRUST}, then \code{P_opt} is optimal \eqn{\delta}.}

\item{ACD_matrix}{a matrix that returns \code{ACD} for different values of a tuning parameter.
If method is \code{DBSCAN}, then the tuning parameter is \eqn{\epsilon}. 
If method is \code{TRUST}, then the tuning parameter is \eqn{\delta}.}
}
              

\note{
Parameters \code{lb,ub} are end points for a range of search for the optimal parameter.
The parameter candidates are calculated in a way such that \eqn{P:=  1.1^x , x \in {lb,lb+0.5,lb+1.0,...,ub }}. Although the default range of search is sufficiently wide, in some cases \code{lb,ub} can be further extended if a warning message is given.

For more discussion on properties of the considered clustering algorithms and DR procedure see Huang et al. (2016).
}
            
\references{

Ciampi, A., Appice, A., and Malerba, D. (2010). Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.

Ester, M., Kriegel, H. P., Sander, J., and Xu, X. (1996). A density-based algorithm for discovering clusters in large spatial databases with noise. In \emph{Proc. KDD 1996} (Vol. 96, No. 34, pp. 226--231).

Huang, X., Iliev, I., Brenning, A., and Gel, Y. R. (2016). Space-Time Clustering with Stability Probe while Riding Downhill. In \emph{KDD 2016 workshop on Mining and Learning from Time Series (MiLeTS 2016)} \url{http://www-bcf.usc.edu/~liu32/milets16/#papers}.
}
            
\author{
Xin Huang, Yulia R. Gel
}           
            
            
\seealso{
\code{\link{CSlideCluster}}, \code{\link[dbscan]{dbscan}}
}
         
            
\examples{
\dontrun{
## example 1
## use iris data to test DR procedure

data(iris)  
require(clue)  # calculate NMI to compare the clustering result with the ground truth
require(scatterplot3d)

Data <- scale(iris[,-5])
ground_truth_label <- iris[,5]

# perform DR procedure to select optimal eps for DBSCAN 
# and save it in variable eps_opt
eps_opt <- DR(t(Data), method="DBSCAN", minPts = 5)$P_opt   

# apply DBSCAN with the optimal eps on iris data 
# and save the clustering result in variable res
res <- dbscan(Data, eps = eps_opt, minPts =5)$cluster  

# calculate NMI to compare the clustering result with the ground truth label
clue::cl_agreement(as.cl_partition(ground_truth_label),
            as.cl_partition(as.numeric(res)), method = "NMI") 
# visualize the clustering result and compare it with the ground truth result
# 3D visualization of clustering result using variables Sepal.Width, Sepal.Length, 
# and Petal.Length
scatterplot3d(Data[,-4],color = res)
# 3D visualization of ground truth result using variables Sepal.Width, Sepal.Length,
# and Petal.Length
scatterplot3d(Data[,-4],color = as.numeric(ground_truth_label))


## example 2
## use synthetic time series data to test DR procedure

require(funtimes)
require(clue) 
require(zoo)

# simulate 16 time series for 4 clusters, each cluster contains 4 time series
set.seed(114) 
samp_Ind <- sample(12,replace=F)
time_points <- 30
X <- matrix(0,nrow=time_points,ncol = 12)
cluster1 <- sapply(1:4,function(x) arima.sim(list(order=c(1,0,0),ar=c(0.2)),
                                             n=time_points,mean=0,sd=1))
cluster2 <- sapply(1:4,function(x) arima.sim(list(order=c(2,0,0),ar=c(0.1,-0.2)),
                                             n=time_points,mean=2,sd=1))
cluster3 <- sapply(1:4,function(x) arima.sim(list(order=c(1,0,1),ar=c(0.3),ma=c(0.1)),
                                             n=time_points,mean=6,sd=1))

X[,samp_Ind[1:4]] <- t(round(cluster1,4))
X[,samp_Ind[5:8]] <- t(round(cluster2,4))
X[,samp_Ind[9:12]] <- t(round(cluster3,4))


# create ground truth label of the synthetic data
ground_truth_label = matrix(1,nrow=12,ncol=1) 
for(k in 1:3){
  ground_truth_label[samp_Ind[(4*k-4+1):(4*k)]] = k
}

# perform DR procedure to select optimal delta for TRUST
# and save it in variable delta_opt
delta_opt <- DR(X,method="TRUST")$P_opt 

# apply TRUST with the optimal delta on the synthetic data 
# and save the clustering result in variable res
res <- CSlideCluster(X,Delta=delta_opt ,Theta=0.9)  

# calculate NMI to compare the clustering result with the ground truth label
clue::cl_agreement(as.cl_partition(as.numeric(ground_truth_label)),
             as.cl_partition(as.numeric(res)),method = "NMI")

# visualize the clustering result and compare it with the ground truth result
# visualization of the clustering result obtained by TRUST
plot.zoo(X, type = "l",plot.type = "single",col = res, xlab = "Time Index", ylab ="")
# visualization of the ground truth result 
plot.zoo(X, type = "l",plot.type = "single",col = ground_truth_label,
         xlab = "Time Index", ylab ="")

             
}
}

\keyword{ts}
\keyword{trend}
