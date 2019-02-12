#Auxiliary functions related to the TRUST algorithm;
#were in the main body of the package until v.5.


# ##### Neighborhood of Time Series #####
#
# This is an auxiliary function to identify which time series in \code{Bv} 
# are \eqn{E^\theta_\delta}-neighbors of \code{Bu}, based on Definition 2 
# by \insertCite{Ciampi_etal_2010;textual}{funtimes}.
# 
# 
# @param Bu a time series vector for which the neighborhood is investigated.
# @param Bv a time series vector (of the same length as \code{Bu}) or a matrix 
# (time series in columns) containing potential neighbors.
# @param Alpha lower limit of the time series domain.
# @param Beta upper limit of the time series domain.
# @param Delta closeness parameter, a real value in \eqn{[0,1]}.
# @param Theta connectivity parameter, a real value in \eqn{[0,1]}.
# 
# 
# @return A vector of logical values indicating which time series in \code{Bv} 
# are \eqn{E^\theta_\delta}-neighbors of \code{Bu}.
# 
# Bu <- rnorm(10)
# Bv <- rnorm(10)
# Alpha <- min(c(Bu, Bv))
# Beta <- max(c(Bu, Bv))
# CNeighbor(Bu, Bv, Alpha, Beta, Delta = 0.5, Theta = 0.8)

CNeighbor <- function(Bu, Bv, Alpha, Beta, Delta, Theta){
    p <- length(Bu)
    if(is.null(dim(Bv)[2])) {Bv <- matrix(Bv, ncol = 1)}
    colSums(abs(Bu - Bv) / (Beta - Alpha) <= Delta) >= Theta * p
}


# ##### Time Series Cluster Homogeneity #####
# 
# This is an auxiliary function to check homogeneity of time series cluster, based on 
# Definition 4 by \insertCite{Ciampi_etal_2010;textual}{funtimes}.
# 
# 
# @param Bu bucket of time series already included in the cluster.
# @param Bv bucket of time series (neighbors) for potential inclusion in the cluster.
# @inheritParams CNeighbor
# 
# 
# @return A logical value indicating whether time series in \code{Bu} and \code{Bv} 
# form a homogeneous cluster.
# 
# Bu <- rnorm(10)
# Bv <- rnorm(10)
# Alpha <- min(c(Bu, Bv))
# Beta <- max(c(Bu, Bv))
# CHomogeneity(Bu, Bv, Alpha, Beta, Delta = 0.5)

CHomogeneity <- function(Bu, Bv, Alpha, Beta, Delta)
{
    M <- cbind(Bu, Bv)
    med <- apply(M, 2, median)
    MED <- matrix(med, length(med), length(med))
    tmp <- abs(MED - t(MED))/(Beta - Alpha) <= Delta
    all(tmp)
}


# ##### Slide-Level Time Series Cluster Expansion #####
# 
# This is an auxiliary function to expand a slide-level time series cluster, 
# based on \insertCite{Ciampi_etal_2010;textual}{funtimes}.
# 
# 
# @param u a time series vector -- a seed to expand the cluster.
# @param Xuncl a time series vector (of the same length as \code{u}) or a matrix 
# (time series in columns) containing unclustered time series.
# @inheritParams CNeighbor
# 
# 
# @return A vector of logical values indicating which time series in \code{Xuncl} 
# should be included in the slide-level cluster with \code{u}.
# 
# set.seed(123)
# u <- rnorm(10)
# Xuncl <- matrix(rt(50, 5), 10, 5)
# Alpha <- min(cbind(u, Xuncl))
# Beta <- max(cbind(u, Xuncl))
# CExpandSlideCluster(u, Xuncl, Alpha, Beta, Delta = 0.15, Theta = 0.8)

CExpandSlideCluster <- function(u, Xuncl, Alpha, Beta, Delta, Theta)
{
    if(is.null(dim(Xuncl)[2])) {Xuncl <- matrix(Xuncl, ncol = 1)}
    ClusterInclude <- rep(FALSE, dim(Xuncl)[2])
    if(length(ClusterInclude) == 0) {return(ClusterInclude)}
    #Define neighbors of the time series among the unclussified unincluded time series Xuncl
    neighbors <- CNeighbor(u, Xuncl, Alpha, Beta, Delta, Theta)
    if (sum(neighbors) > 0) { #Proceed with all these things only if some neighbors were found
        neighborsID <- which(neighbors)
        #Check the homegeneity of the time series with its all neighbors
        if (CHomogeneity(u, Xuncl[,neighbors], Alpha, Beta, Delta)){ #The seed time series and all its neighbors are homogeneous
            #Add all neighbors to the cluster
            ClusterInclude[neighbors] <- TRUE
            #If all previously unclassified values were included in the cluster, return the result (FINISH).
            if (all(ClusterInclude)) {
                return(ClusterInclude)
            } else {
                #Otherwise expand cluster using each of the neighbors
                for (i in 1:sum(neighbors)){
                    tmp <- CExpandSlideCluster(Xuncl[,neighborsID[i]], 
                                               Xuncl[,!ClusterInclude], 
                                               Alpha, Beta, Delta, Theta)
                    ClusterInclude[!ClusterInclude][tmp] <- TRUE
                }
            }
        } else { #The cluster is not homogeneous
            #Check each neighbor one by one
            for (i in 1:sum(neighbors)){
                if (CHomogeneity(u, Xuncl[,neighborsID[i]], Alpha, Beta, Delta)){
                    ClusterInclude[neighborsID[i]] <- TRUE
                    tmp <- CExpandSlideCluster(Xuncl[,neighborsID[i]], 
                                               Xuncl[,!ClusterInclude], 
                                               Alpha, Beta, Delta, Theta)
                    ClusterInclude[!ClusterInclude][tmp] <- TRUE
                }
            }  
        }
    } 
    return(ClusterInclude)
}


# ##### Window-Level Time Series Cluster Expansion #####
# 
# This is an auxiliary function to expand a window-level time series cluster, 
# based on \insertCite{Ciampi_etal_2010;textual}{funtimes}.
# 
# 
# @param e a vector of logical values identifying which time series among \code{Euncl} 
# were clustered together with \code{e} over at least \code{w*Epsilon} slides 
# within a window; see Definition 7 by \insertCite{Ciampi_etal_2010;textual}{funtimes}. 
# This is a seed for window-level clustering.
# @param Euncl a square matrix identifying the binary window cluster relation 
# for yet unclustered time series.
# 
# 
# @return A vector of logical values indicating which time series in \code{Euncl} 
# should be included in the window-level cluster with \code{e}.
# 
# set.seed(123)
# e <- sample(c(TRUE, FALSE), 5, replace = TRUE)
# Euncl <- matrix(sample(c(TRUE, FALSE), 5, replace = TRUE), 5, 5)
# CExpandWindowCluster(e, Euncl)

CExpandWindowCluster <- function(e, Euncl)
{
    if(is.null(dim(Euncl)[2])) {Euncl <- matrix(Euncl, ncol = 1)}
    ClusterInclude <- rep(FALSE, dim(Euncl)[2])
    if(length(ClusterInclude) == 0) {return(ClusterInclude)}
    if(any(e)) { #Perform it only if we have neighbors
        ClusterInclude[e] <- TRUE
        NeighID <- which(e)
        for (i in 1:sum(e)) {
            if (all(ClusterInclude)) {return(ClusterInclude)}
            tmp <- CExpandWindowCluster(Euncl[!ClusterInclude,NeighID[i]], 
                                        Euncl[!ClusterInclude, !ClusterInclude])
            ClusterInclude[!ClusterInclude][tmp] <- TRUE
        }
    }  
    return(ClusterInclude)
}

