#' @importFrom dbscan dbscan

.para_opt <- function(M) {
    acd <- M[, 2]
    max_acd <- max(acd)
    if (max_acd == 0) {
        warning("All ACDs are 0. The range of controlling parameters might be unsuitable.")
        return(NA)
    }
    
    # Find the last index of the maximum ACD value
    max_indices <- which(acd == max_acd)
    max_idx <- max_indices[length(max_indices)]
    
    # Check if there is enough room after the max ACD to find a local minimum.
    # The logic looks for a minimum in a window of size 3, so we need at least 3 points after max_idx+3.
    if (max_idx > (length(acd) - 6)) {
        warning("Max ACD is found near the upper bound. Consider increasing the 'ub' parameter.")
        return(NA)
    }
    
    # Search for the first local minimum after the max ACD.
    # A point is a local minimum if it's smaller than or equal to its 3 neighbors to the left,
    # and strictly smaller than its 3 neighbors to the right.
    for (i in (max_idx + 3):(length(acd) - 3)) {
        is_local_min <- all(acd[i] <= acd[(i - 3):(i - 1)]) && 
                        all(acd[i] < acd[(i + 1):(i + 3)])
        if (is_local_min) {
            return(M[i, 1])
        }
    }
    
    warning("Could not find a suitable local minimum for the optimal parameter.")
    return(NA)
}

.calculate_acd <- function(param, data, B_samples, n_nodes, cluster_func) {
    diffs <- replicate(B_samples, {
        s_index <- sample(ncol(data), n_nodes, replace = FALSE)
        cl1 <- max(cluster_func(data[, s_index, drop = FALSE], param))
        cl2 <- max(cluster_func(data[, -s_index, drop = FALSE], param))
        abs(cl1 - cl2)
    })
    mean(diffs)
}

#' Downhill Riding (DR) Procedure
#' 
#' Downhill riding procedure for selecting optimal tuning parameters in clustering 
#' algorithms, using an (in)stability probe.
#' 
#' @details Parameters \code{lb,ub} are endpoints for the search for the 
#' optimal parameter. The parameter candidates are calculated in a way such that 
#' \eqn{P:=  1.1^x , x \in {lb,lb+0.5,lb+1.0,...,ub}}. 
#' Although the default range of search is sufficiently wide, in some cases 
#' \code{lb,ub} can be further extended if a warning message is given.
#' 
#' For more discussion on properties of the considered clustering algorithms and the 
#' DR procedure see \insertCite{Huang_etal_2016;textual}{funtimes} 
#' and \insertCite{Huang_etal_2018_riding;textual}{funtimes}.
#' 
#' 
#' @param X an \eqn{n	imes k} matrix where columns are \eqn{k} objects to be clustered, 
#' and each object contains n observations (objects could be a set of time series).
#' @param method the clustering method to be used -- currently either 
#' \dQuote{TRUST} \insertCite{Ciampi_etal_2010}{funtimes} 
#' or \dQuote{DBSCAN} \insertCite{Ester_etal_1996}{funtimes}. If the method is \code{DBSCAN}, 
#' then set \code{MinPts} and optimal \eqn{\epsilon} is selected using DR. 
#' If the method is \code{TRUST}, then set \code{theta}, and optimal \eqn{\delta} 
#' is selected using DR.
#' @param minPts the minimum number of samples in an \eqn{\epsilon}-neighborhood of 
#' a point to be considered as a core point. The \code{minPts} is to be used only 
#' with the \code{DBSCAN} method. The default value is 3.
#' @param theta connectivity parameter \eqn{	heta \in (0,1)}, which is to be used 
#' only with the \code{TRUST} method. The default value is 0.9.
#' @param B number of random splits in calculating the 
#' Average Cluster Deviation (ACD). The default value is 500.
#' @param lb,ub endpoints for a range of search for the optimal parameter.
#' 
#' 
#' @return A list containing the following components:
#' \item{P_opt}{the value of the optimal parameter. If the method is \code{DBSCAN}, then 
#' \code{P_opt} is optimal \eqn{\epsilon}. If the method is \code{TRUST}, 
#' then \code{P_opt} is optimal \eqn{\delta}.}
#' \item{ACD_matrix}{a matrix that returns \code{ACD} for different values of a 
#' tuning parameter.
#' If the method is \code{DBSCAN}, then the tuning parameter is \eqn{\epsilon}. 
#' If the method is \code{TRUST}, then the tuning parameter is \eqn{\delta}.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{BICC}}, \code{\link[dbscan]{dbscan}}
#' 
#' @keywords ts trend
#' 
#' @author Xin Huang, Yulia R. Gel
#' 
#' @export
#' @examples
#' \dontrun{
#' ## example 1
#' ## use iris data to test DR procedure
#' 
#' data(iris)  
#' require(clue)  # calculate NMI to compare the clustering result with the ground truth
#' require(scatterplot3d)
#' 
#' Data <- scale(iris[,-5])
#' ground_truth_label <- iris[,5]
#' 
#' # perform DR procedure to select optimal eps for DBSCAN 
#' # and save it in variable eps_opt
#' eps_opt <- DR(t(Data), method="DBSCAN", minPts = 5)$P_opt   
#' 
#' # apply DBSCAN with the optimal eps on iris data 
#' # and save the clustering result in variable res
#' res <- dbscan(Data, eps = eps_opt, minPts =5)$cluster  
#' 
#' # calculate NMI to compare the clustering result with the ground truth label
#' clue::cl_agreement(as.cl_partition(as.numeric(ground_truth_label)),
#'                    as.cl_partition(as.numeric(res)), method = "NMI") 
#' # visualize the clustering result and compare it with the ground truth result
#' # 3D visualization of clustering result using variables Sepal.Width, Sepal.Length, 
#' # and Petal.Length
#' scatterplot3d(Data[,-4],color = res)
#' # 3D visualization of ground truth result using variables Sepal.Width, Sepal.Length,
#' # and Petal.Length
#' scatterplot3d(Data[,-4],color = as.numeric(ground_truth_label))
#' 
#' 
#' ## example 2
#' ## use synthetic time series data to test DR procedure
#' 
#' require(funtimes)
#' require(clue) 
#' require(zoo)
#' 
#' # simulate 16 time series for 4 clusters, each cluster contains 4 time series
#' set.seed(114) 
#' samp_Ind <- sample(12,replace=F)
#' time_points <- 30
#' X <- matrix(0,nrow=time_points,ncol = 12)
#' cluster1 <- sapply(1:4,function(x) arima.sim(list(order = c(1, 0, 0), ar = c(0.2)),
#'                                              n = time_points, mean = 0, sd = 1))
#' cluster2 <- sapply(1:4,function(x) arima.sim(list(order = c(2 ,0, 0), ar = c(0.1, -0.2)),
#'                                              n = time_points, mean = 2, sd = 1))
#' cluster3 <- sapply(1:4,function(x) arima.sim(list(order = c(1, 0, 1), ar = c(0.3), ma = c(0.1)),
#'                                              n = time_points, mean = 6, sd = 1))
#' 
#' X[,samp_Ind[1:4]] <- t(round(cluster1, 4))
#' X[,samp_Ind[5:8]] <- t(round(cluster2, 4))
#' X[,samp_Ind[9:12]] <- t(round(cluster3, 4))
#' 
#' 
#' # create ground truth label of the synthetic data
#' ground_truth_label = matrix(1, nrow = 12, ncol = 1) 
#' for(k in 1:3){
#'     ground_truth_label[samp_Ind[(4*k - 4 + 1):(4*k)]] = k
#' }
#' 
#' # perform DR procedure to select optimal delta for TRUST
#' # and save it in variable delta_opt
#' delta_opt <- DR(X, method = "TRUST")$P_opt 
#' 
#' # apply TRUST with the optimal delta on the synthetic data 
#' # and save the clustering result in variable res
#' res <- CSlideCluster(X, Delta = delta_opt, Theta = 0.9)  
#' 
#' # calculate NMI to compare the clustering result with the ground truth label
#' clue::cl_agreement(as.cl_partition(as.numeric(ground_truth_label)),
#'                    as.cl_partition(as.numeric(res)), method = "NMI")
#' 
#' # visualize the clustering result and compare it with the ground truth result
#' # visualization of the clustering result obtained by TRUST
#' plot.zoo(X, type = "l", plot.type = "single", col = res, xlab = "Time index", ylab = "")
#' # visualization of the ground truth result 
#' plot.zoo(X, type = "l", plot.type = "single", col = ground_truth_label,
#'          xlab = "Time index", ylab = "")
#' }
#' 
DR <- function(X, method, minPts = 3, theta = 0.9, B = 500, lb = -30, ub = 10)
{
    control_para <- 1.1^seq(lb, ub, by = 0.5)
    Nnodes <- floor(ncol(X) / 2)

    message(sprintf("Total # of controlling parameters: %d", length(control_para)))
    
    if (method == "DBSCAN") {
        if (missing(minPts)) {
            stop("Missing parameter minPts for DBSCAN")
        }
        dbscan_fn <- function(data, eps) dbscan::dbscan(t(data), eps = eps, minPts = minPts)$cluster
        acd_values <- sapply(seq_along(control_para), function(i) {
            message(sprintf("Processing parameter %d of %d", i, length(control_para)))
            .calculate_acd(control_para[i], X, B, Nnodes, dbscan_fn)
        })
    } else if (method == "TRUST") {
        trust_fn <- function(data, delta) CSlideCluster(data, Delta = delta, Theta = theta)
        acd_values <- sapply(seq_along(control_para), function(i) {
            message(sprintf("Processing parameter %d of %d", i, length(control_para)))
            .calculate_acd(control_para[i], X, B, Nnodes, trust_fn)
        })
    } else {
        stop("Method must be 'DBSCAN' or 'TRUST'")
    }

    ACD_matrix <- data.frame("Controlling.Parameter" = control_para, "ACD" = acd_values)
    
    p_opt <- .para_opt(ACD_matrix)
    
    out <- list(P_opt = p_opt, ACD_matrix = ACD_matrix)
    return(out)
}
