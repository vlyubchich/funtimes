#' Window-Level Time Series Clustering
#' 
#' Cluster time series at a window level, 
#' based on Algorithm 2 of \insertCite{Ciampi_etal_2010;textual}{funtimes}.
#' 
#' @details This is the upper-level function for time series clustering. It exploits 
#' the function \code{\link{CSlideCluster}} to cluster time series within each slide 
#' based on closeness and homogeneity measures. Then, it uses slide-level cluster 
#' assignments to cluster time series within each window.
#' 
#' The total length of time series (number of levels, i.e., \code{nrow(X)}) 
#' should be divisible by \code{p}.
#' 
#' 
#' @inheritParams CSlideCluster
#' @param Alpha lower limit of the time series domain, 
#' passed to \code{\link{CSlideCluster}}.
#' @param Beta upper limit of the time series domain, passed to \code{\link{CSlideCluster}}.
#' @param Delta closeness parameter, passed to \code{\link{CSlideCluster}}.
#' @param Theta connectivity parameter, passed to \code{\link{CSlideCluster}}.
#' @param p number of layers (time series observations) in each slide.
#' @param w number of slides in each window.
#' @param s step to shift a window, calculated in number of slides. The recommended 
#' values are 1 (overlapping windows) or equal to \code{w} (non-overlapping windows).
#' @param Epsilon a real value in \eqn{[0,1]} used to identify each pair of time series 
#' that are clustered together over at least \code{w*Epsilon} slides within a window;  
#' see Definition 7 by \insertCite{Ciampi_etal_2010;textual}{funtimes}. Default is 1.
#' 
#' 
#' @return A vector (if \code{X} contains only one window) or matrix with cluster 
#' labels for each time series (columns) and window (rows). 
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{CSlideCluster}}, \code{\link{CWindowCluster}}, 
#' and \code{\link{BICC}}
#' 
#' @keywords cluster ts trend
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' #For example, weekly data come in slides of 4 weeks
#' p <- 4 #number of layers in each slide (data come in a slide)
#'     
#' #We want to analyze the trend clusters within a window of 1 year
#' w <- 13 #number of slides in each window
#' s <- w  #step to shift a window
#' 
#' #Simulate 26 autoregressive time series with two years of weekly data (52*2 weeks), 
#' #with a 'burn-in' period of 300.
#' N <- 26
#' T <- 2*p*w
#'     
#' set.seed(123) 
#' phi <- c(0.5) #parameter of autoregression
#' X <- sapply(1:N, function(x) arima.sim(n = T + 300, 
#'      list(order = c(length(phi), 0, 0), ar = phi)))[301:(T + 300),]
#' colnames(X) <- paste("TS", c(1:dim(X)[2]), sep = "")
#'  
#' tmp <- CWindowCluster(X, Delta = NULL, Theta = 0.8, p = p, w = w, s = s, Epsilon = 1)
#' 
#' #Time series were simulated with the same parameters, but based on the clustering parameters,
#' #not all time series join the same cluster. We can plot the main cluster for each window, and 
#' #time series out of the cluster:
#' par(mfrow = c(2, 2))
#' ts.plot(X[c(1:(p*w)), tmp[1,] == 1], ylim = c(-4, 4), 
#'         main = "Time series cluster 1 in window 1")
#' ts.plot(X[c(1:(p*w)), tmp[1,] != 1], ylim = c(-4, 4), 
#'         main = "The rest of the time series in window 1")
#' ts.plot(X[c(1:(p*w)) + s*p, tmp[2,] == 1], ylim = c(-4, 4), 
#'         main = "Time series cluster 1 in window 2")
#' ts.plot(X[c(1:(p*w)) + s*p, tmp[2,] != 1], ylim = c(-4, 4), 
#'         main = "The rest of the time series in window 2")
#' 
CWindowCluster <- function(X, Alpha = NULL, Beta = NULL, Delta = NULL, Theta = 0.8, 
                           p, w, s, Epsilon = 1)
{
    T <- dim(X)[1]
    N <- dim(X)[2]
    nWindows <- length(seq(from = p*w, to = T, by = s*p))
    Clusters <- array(NA, dim = c(w, N, nWindows))
    ClustersWithinWindow <- array(NA, dim = c(nWindows, N))
    #Separate data into windows
    for (nw in 1:nWindows) {
        #Apply CSlideCluster to each slide withing the window
        for (sl in 1:w) {
            Clusters[sl,,nw] <- CSlideCluster(X[c(((nw - 1)*w*p + 1 + (sl - 1)*p):((nw - 1)*w*p + (sl)*p)),], 
                                              Alpha = Alpha, Beta = Beta, Delta = Delta, Theta = Theta)
        }
        #Apply clusterization to CSlideCluster within window
        #count how many times each node was clustered together with each other node, frow w times
        E <- (sapply(1:N, function(n) colSums(Clusters[,n,nw] == Clusters[,,nw])) >= Epsilon * w)
        TSclusters <- rep(NA, N)
        currentTS <- currentCluster <- 1
        TSclusters[currentTS] <- currentCluster
        Unclassified <- is.na(TSclusters)
        while (any(Unclassified)) {
            Include <- CExpandWindowCluster(E[Unclassified, currentTS], 
                                            E[Unclassified, Unclassified])
            #Assign cluster number to the time series that were just grouped
            TSclusters[Unclassified][Include] <- currentCluster
            #Next cluster number will be:
            currentCluster <- currentCluster + 1
            #Still unclassified time series
            Unclassified <- is.na(TSclusters)
            #Select the next time series to start the cluster with, and immediately
            #assign it to the new cluster
            currentTS <- which(Unclassified)[1]
            TSclusters[currentTS] <- currentCluster
            Unclassified <- is.na(TSclusters)
        }
        ClustersWithinWindow[nw,] <- TSclusters
    }
    return(ClustersWithinWindow)
}
