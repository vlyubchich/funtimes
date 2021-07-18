#' Topological Clustering using Betti Numbers (TopoCBN)
#' 
#' The function performs unsupervised clustering of multivariate data based on topological data
#' analysis (TDA). The objective is to partition data into non-overlapping clusters, where the 
#' definition of a cluster falls under #' a general framework of density based clustering, 
#' e.g., DBSCAN, OPTICS etc. That is, intuitively the cluster is to be a subset of points which
#' is path-connected, i.e. any point in the subset can be reached from any other one through 
#' a path consisting of points (also belonging to the subset); furthermore, the consecutive 
#' points on the path are close enough and their local neighborhoods are similar in 
#' shape \insertCite{islambekov2019unsupervised}{funtimes}. To compare shapes, TopoCBN builds 
#' a Vietoris-Rips (VR) filtration upon such neighborhoods around each point and compute 
#' topological summaries in the form of the Betti sequences using persistent homology. The closer
#' the Betti sequences to one another for a pair of close-by points, the more likely similar the 
#' shapes of their neighborhoods. Thus, when identifying clusters, TopoCBN utilizes both the
#' distance function and local geometric information around the points. Note that accounting 
#' for shape similarity can be viewed as an extension of conventional clustering properties in 
#' the density-based clustering framework.
#' 
#' @param data A point cloud given as an N by d matrix, where N = number of points, d = dimension of Euclidean space
#' OR an N by N matrix of pairwise distances.
#' @param nKNN Number of k nearest neighbors to take around each point.
#' @param filt_len Filtration length (also length of Betti sequences). Default is 25.
#' @param dist_matrix  Is set to FALSE by default, if data is a point cloud. dist_matrix=TRUE if data
#' is a matrix of pairwise distances.
#' 
#' @return A list with the following components:
#' \item{assignments}{Cluster labels (vector of length N)}
#' \item{nClust}{Number of clusters}
#' \item{cSize}{Cluster sizes (vector of length nClust)}
#' 
#' @import igraph TDA randomcoloR graphics
#' @importFrom FNN knn.index knn.dist
#' @references 
#' \insertAllCited{}
#' 
#' @seealso \code{\link{cumsumCPA_test}} for change point detection in time series via a linear
#' regression with temporally correlated errors
#' 
#' @keywords cluster topology Betti
#' 
#' @author Palina Niamkova, Umar Islambekov, Yulia R. Gel
#' @export 
#' @examples
#' \dontrun{
#' #Example 1:
#' #Let's import dataset with today's Covid-19 parameters per each state:
#' data<-covid19us::get_states_current()
#' #For this example we will keep data for positive cases and deaths today:
#' data<-data[c(3,9)]
#' #We also need to replace NA values to integer 0:  
#' data[is.na(data)] = 0
#' 
#' #Now we can run CBN:
#' result <- TopoCBN(data,nKNN=12) # can also try with filt_len=50,75,100
#' 
#' #We can obtain the same results using matrix of pairwise distances:
#' dMatrix <- as.matrix(dist(data))
#' result <- TopoCBN(dMatrix,nKNN=12,dist_matrix = TRUE)
#' 
#' #Let's plot the results:
#' set.seed(365)
#' distinct_clrs=randomcoloR::distinctColorPalette(result$nClust)
#' clrs<-distinct_clrs[result$assignments] # distinct colors for clusters
#' plot(data,col=clrs,pch=20,xlab='x',ylab='y',main = 'TopoCBN') 
#' print(result)
#' 
#' #We can see that CBN function identified 6 clusters within our dataset.
#'
#' #Example 2:
#' #Let's import dataset with air quality level in  Californian metropolitan areas. The three
#' #columns of the dataset contains indicator of air quality (the lower the better), value 
#' #added of companies (in thousands of dollars).
#' 
#' data<-as.matrix(Ecdat::Airq[1:3])
#' 
#' #Now apply TopoCBN function to the air quality data:
#' result <- TopoCBN(data,nKNN=12) # can also try with filt_len=50,75,100
#'
#' #The same results can be obtained using matrix of pairwise distances:
#' dMatrix <- as.matrix(dist(data))
#' result <- TopoCBN(dMatrix,nKNN=12,dist_matrix = TRUE)
#' 
#' #Plot the results:
#' set.seed(365)
#' distinct_clrs=randomcoloR::distinctColorPalette(result$nClust)
#' clrs<-distinct_clrs[result$assignments] # distinct colors for clusters
#' plot(data,col=clrs,pch=20,xlab='x',ylab='y',main = 'TopoCBN') 
#' print(result)
#' 
#' #We see that TopoCBN identified 4 clusters within our dataset of the sizes 
#' #1,3,5, and 21. These results suggest that companies with added values under $5,000 may 
#' #have any value of air pollution. However, companies with higher added values (>$5,000)
#' #correspond to the dramatically increased (deteriorated) levels of air pollution.
#' }
#' 
TopoCBN = function(data,nKNN,filt_len=25,dist_matrix=FALSE){
  
  # extract betti sequence from PD
  extract_betti=function(PD){
    b<-numeric(length = filt_len)
    for (k in 1:filt_len)
      b[k]<-sum((scale_seq[k]>=PD[,1])&(scale_seq[k]<PD[,2]))
    b
  }
  #main body#
  
  N <- nrow(data) # number of data points
  if (dist_matrix){
    ind <- t(apply(data,1,order))[,1:nKNN]
    maxscale <- max(t(apply(data,1,sort))[,nKNN])
  } else{
    ind <- cbind(1:N,FNN::knn.index(data,k=nKNN-1)) # indices of k nearest neighbors
    maxscale<-max(FNN::knn.dist(data,k=nKNN-1)[,nKNN-1]) # maxscale
  }
  
  scale_seq<-seq(0,maxscale,length.out = filt_len) # increasing sequence of scale values
  
  # Compute betti sequences
  print("Extracting Betti sequences (1/4)")
  betti_0<-betti_1<-matrix(nrow = N,ncol = filt_len) 
  for (i in 1:N){
    
    # construct Rips filtration built on neighborhood around point i
    if (dist_matrix){
      ripsFltr<-TDA::ripsFiltration(data[ind[i,],ind[i,]],maxdimension = 0,maxscale = maxscale,dist = 'arbitrary') 
    } else{
      ripsFltr<-TDA::ripsFiltration(data[ind[i,],],maxdimension = 0,maxscale = maxscale,dist = 'euclidean') 
    }
    # compute PD of Rips filtration
    PD<-TDA::filtrationDiag(filtration = ripsFltr, maxdimension = 1)$diagram
    # replace death=Inf with death=maxscale
    PD[PD[,3]==Inf,3]=maxscale
    
    # extract betti-0 and betti-1 sequences from PD
    betti_0[i,]=extract_betti(PD[PD[,1]==0,2:3,drop=F])
    betti_1[i,]=extract_betti(PD[PD[,1]==1,2:3,drop=F])
    
  }
  
  ##################################################
  # Computing relative change in betti sequences
  ##################################################
  print("Computing relative change in Betti sequences (2/4)")
  delta_betti_0<-delta_betti_1<-matrix(nrow = N,ncol = nKNN)
  for (i in 1:N){
    
    bi_bj<-as.matrix(dist(betti_0[ind[i,],]))[1,]
    denom<-norm(as.matrix(betti_0[i,]),type = 'f')
    delta_betti_0[i,]=bi_bj/denom
    
    bi_bj<-as.matrix(dist(betti_1[ind[i,],]))[1,]
    denom<-norm(as.matrix(betti_1[i,]),type = 'f')
    delta_betti_1[i,]=bi_bj/denom 
  }
  
  ##################################################
  # computing cutoff thresholds
  ##################################################
  
  bp_0<-boxplot(as.vector(delta_betti_0),plot = F)
  bp_1<-boxplot(as.vector(delta_betti_1),plot = F)
  
  cutoff_0<-bp_0$stats[5,] # cutoff threshold for betti-0. Can also be set manually
  cutoff_1<-bp_1$stats[5,] # cutoff threshold for betti-1. Can also be set manually
  
  ############################
  # Forming adjacency matrix
  ############################
  print('Forming adjacency matrix (3/4)')
  A<-matrix(0,ncol = N,nrow = N) 
  for (i in 1:N)
  {
    index_0<-ind[i,which(delta_betti_0[i,]<=cutoff_0)]
    index_1<-ind[i,which(delta_betti_1[i,]<=cutoff_1)]
    
    A[i,intersect(index_0,index_1)]=1
    
  }
  ############################
  # Performing clustering
  ############################
  print('Performing clustering (4/4)')
  g<-igraph::graph_from_adjacency_matrix(A,mode='directed') # form graph g from adj. matrix A
  clstrs<-clusters(g,mode = 'strong') # define clusters as strongly connected compents of graph g
  
  return(list(assignments=clstrs$membership,nClust=clstrs$no,cSize=clstrs$csize))
}


