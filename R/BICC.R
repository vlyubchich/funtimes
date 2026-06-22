#' BIC-Based Spatio-Temporal Clustering
#'
#' Apply the algorithm of unsupervised spatio-temporal clustering, TRUST 
#' \insertCite{Ciampi_etal_2010}{funtimes}, with automatic selection of its 
#' tuning parameters \code{Delta} and \code{Epsilon} based on Bayesian 
#' information criterion, BIC \insertCite{Schaeffer_etal_2016_trust}{funtimes}.
#'
#' @details This is the upper-level function for time series clustering. 
#' It exploits the functions \code{\link{CWindowCluster}} and 
#' \code{\link{CSlideCluster}} to cluster time series based on closeness and 
#' homogeneity measures. Clustering is performed multiple times with a range 
#' of equidistant values for the parameters \code{Delta} and \code{Epsilon}, 
#' then optimal parameters \code{Delta} and \code{Epsilon} along with the 
#' corresponding clustering results are shown 
#' \insertCite{@see @Schaeffer_etal_2016_trust, for more details}{funtimes}.
#' 
#' The total length of time series (number of levels, i.e., \code{nrow(X)}) 
#' should be divisible by \code{p}.
#' 
#' 
#' @inheritParams CWindowCluster
#' 
#' 
#' @return A list with the following elements:
#' \item{delta.opt}{optimal value for the clustering parameter \code{Delta}.}
#' \item{epsilon.opt}{optimal value for the clustering parameter \code{Epsilon}.}
#' \item{clusters}{vector of length \code{ncol(X)} with cluster labels.}
#' \item{IC}{values of the information criterion (BIC) for each considered 
#' combination of \code{Delta} (rows) and \code{Epsilon} (columns).}
#' \item{delta.all}{vector of considered values for \code{Delta}.}
#' \item{epsilon.all}{vector of considered values for \code{Epsilon}.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{CSlideCluster}}, \code{\link{CWindowCluster}}, \code{\link{purity}}
#' 
#' @keywords cluster ts trend
#' 
#' @author Ethan Schaeffer, Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' # Fix seed for reproducible simulations:
#' set.seed(1)
#' 
#' ##### Example 1
#' # Similar to Schaeffer et al. (2016), simulate 3 years of monthly data 
#' #for 10 locations and apply clustering:
#' # 1.1 Simulation
#' T <- 36 #total months
#' N <- 10 #locations
#' phi <- c(0.5) #parameter of autoregression
#' burn <- 300 #burn-in period for simulations
#' X <- sapply(1:N, function(x) 
#'     arima.sim(n = T + burn, 
#'               list(order = c(length(phi), 0, 0), ar = phi)))[(burn + 1):(T + burn),]
#' colnames(X) <- paste("TS", c(1:dim(X)[2]), sep = "")
#' 
#' # 1.2 Clustering
#' # Assume that information arrives in year-long slides or data chunks
#' p <- 12 #number of time layers (months) in a slide
#' # Let the upper level of clustering (window) be the whole period of 3 years, so
#' w <- 3 #number of slides in a window
#' s <- w #step to shift a window, but it does not matter much here as we have only one window of data
#' tmp <- BICC(X, p = p, w = w, s = s)
#' 
#' # 1.3 Evaluate clustering
#' # In these simulations, it is known that all time series belong to one class,
#' #since they were all simulated the same way:
#' classes <- rep(1, 10)
#' # Use the information on the classes to calculate clustering purity:
#' purity(classes, tmp$clusters[1,])
#' 
#' ##### Example 2
#' # 2.1 Modify time series and update classes accordingly:
#' # Add a mean shift to a half of the time series:
#' X2 <- X
#' X2[, 1:(N/2)] <- X2[, 1:(N/2)] + 3
#' classes2 <- rep(c(1, 2), each = N/2)
#' 
#' # 2.2 Re-apply clustering procedure and evaluate clustering purity:
#' tmp2 <- BICC(X2, p = p, w = w, s = s)
#' tmp2$clusters
#' purity(classes2, tmp2$clusters[1,])
#' 
BICC <- function(X, Alpha = NULL, Beta = NULL, Theta = 0.8, p, w, s) 
{
    dims <- dim(X)
    N <- dims[2] #number of time series
    n <- dims[1] #length of time series
    
    IQRx <- median(apply(X, 2, IQR))
    deltas <- seq(IQRx/50, IQRx/2, IQRx/10)
    epsilons <- seq(1.0/w, 1.0, 1.0/w)
    
    IC <- array(NA, dim = c(length(deltas), length(epsilons)))
    
    for (i in 1:length(deltas)) { 
        for (j in 1:length(epsilons)) {
            clusters_tmp <- CWindowCluster(X, Delta = deltas[i], Epsilon = epsilons[j], 
                                           Alpha = Alpha, Beta = Beta, Theta = Theta, 
                                           p = p, w = w, s = s)
            k <- max(clusters_tmp) 
            
            # Calculate residuals for each cluster
            residuals_list <- lapply(1:k, function(v) {
                cluster_series <- X[, clusters_tmp == v, drop = FALSE]
                cluster_mean <- rowMeans(cluster_series)
                cluster_series - cluster_mean
            })
            
            # Biased variance of all residuals
            all_residuals <- unlist(residuals_list)
            residuals_var <- var(all_residuals) * (n * N - 1) / (n * N)
            
            IC[i, j] <- N * log(residuals_var) + k * log(N) #BIC
        }
    }
    
    IC[is.infinite(IC)] <- NA
    
    # Find optimal parameters.
    # In case of ties for the minimum IC, which.min returns the first one.
    # We are following that behavior by taking the first row of the index matrix.
    min_ic_idx <- which(IC == min(IC, na.rm = TRUE), arr.ind = TRUE)
    opt_idx <- min_ic_idx[1, ]
    delta_opt <- deltas[opt_idx[1]]
    epsilon_opt <- epsilons[opt_idx[2]]
    
    output <- CWindowCluster(X, Delta = delta_opt, Epsilon = epsilon_opt, 
                             Alpha = Alpha, Beta = Beta, Theta = Theta, 
                             p = p, w = w, s = s)
    
    dimnames(output) <- list(paste("Window", c(1:nrow(output)), sep = "_"), colnames(X))
    
    return(list(delta.opt = delta_opt, 
                epsilon.opt = epsilon_opt, 
                clusters = output, 
                IC = IC, 
                delta.all = deltas, 
                epsilon.all = epsilons))
}
