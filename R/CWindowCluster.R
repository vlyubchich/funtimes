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
#' @param Alpha lower limit of the time-series domain, 
#' passed to \code{\link{CSlideCluster}}.
#' @param Beta upper limit of the time-series domain passed to \code{\link{CSlideCluster}}.
#' @param Delta closeness parameter passed to \code{\link{CSlideCluster}}.
#' @param Theta connectivity parameter passed to \code{\link{CSlideCluster}}.
#' @param p number of layers (time-series observations) in each slide.
#' @param w number of slides in each window.
#' @param s step to shift a window, calculated in the number of slides. The recommended 
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
    X <- as.matrix(X)
    if (!is.numeric(X))
        stop("X must be numeric.")
    if (any(is.na(X)))
        stop("X contains missing values.")
    if (ncol(X) < 2)
        stop("X must contain at least two time series (columns).")

    p <- as.integer(p)
    w <- as.integer(w)
    s <- as.integer(s)
    if (any(is.na(c(p, w, s))) || any(c(p, w, s) < 1))
        stop("p, w, and s must be positive integers.")
    if (!is.numeric(Epsilon) || length(Epsilon) != 1L || is.na(Epsilon) || Epsilon < 0 || Epsilon > 1)
        stop("Epsilon must be a single numeric value in [0, 1].")

    T_len <- nrow(X)
    N <- ncol(X)
    if (T_len %% p != 0)
        stop("nrow(X) must be divisible by p.")
    if (p * w > T_len)
        stop("p * w must be less than or equal to nrow(X).")
    
    n_windows <- length(seq(from = p * w, to = T_len, by = s * p))
    if (n_windows < 1)
        stop("No windows can be formed with the provided p, w, and s values.")
    
    slide_clusters <- array(NA, dim = c(w, N, n_windows))
    window_clusters <- array(NA, dim = c(n_windows, N))
    
    #Separate data into windows
    for (nw in 1:n_windows) {
        # Apply CSlideCluster to each slide within the window
        for (sl in 1:w) {
            slide_start_idx <- (nw - 1) * s * p + (sl - 1) * p + 1
            slide_end_idx <- (nw - 1) * s * p + sl * p
            slide_data <- X[slide_start_idx:slide_end_idx, , drop = FALSE]
            slide_clusters[sl, , nw] <- CSlideCluster(slide_data, 
                                                      Alpha = Alpha, Beta = Beta, Delta = Delta, Theta = Theta)
        }
        
        # E is a similarity matrix: E[i,j] is TRUE if TS i and j are often in the same slide-cluster.
        E <- (sapply(1:N, function(n) colSums(slide_clusters[, n, nw] == slide_clusters[, , nw])) >= Epsilon * w)
        
        # Cluster based on the similarity matrix E
        ts_clusters_in_window <- rep(NA, N)
        cluster_label <- 1
        unclassified_indices <- 1:N
        
        while (length(unclassified_indices) > 0) {
            seed_ts_idx <- unclassified_indices[1]
            ts_clusters_in_window[seed_ts_idx] <- cluster_label
            
            unclassified_subset_indices <- unclassified_indices[-1]
            
            if (length(unclassified_subset_indices) > 0) {
                seed_similarities <- E[unclassified_subset_indices, seed_ts_idx]
                subset_similarity_matrix <- E[unclassified_subset_indices, unclassified_subset_indices, drop = FALSE]
                
                series_to_include <- CExpandWindowCluster(seed_similarities, subset_similarity_matrix)
                
                ts_clusters_in_window[unclassified_subset_indices[series_to_include]] <- cluster_label
            }
            
            cluster_label <- cluster_label + 1
            unclassified_indices <- which(is.na(ts_clusters_in_window))
        }
        window_clusters[nw, ] <- ts_clusters_in_window
    }
    return(window_clusters)
}
