#' Time Series Clustering based on Trend Synchronism
#' 
#' Cluster time series with a common parametric trend using the 
#' \code{\link{sync.test}} function 
#' \insertCite{Lyubchich_Gel_2016_synchronism,Ghahari_etal_2017_MBDCE}{funtimes}.
#' 
#' @details The \code{sync.cluster} function recursively clusters time series having 
#' a pre-specified common parametric trend until there are no time series left. 
#' Starting with the given \eqn{N} time series, the \code{\link{sync.test}} function 
#' is used to test for a common trend. If null hypothesis of common trend is not 
#' rejected by \code{\link{sync.test}}, the time series are grouped together 
#' (i.e., assigned to a cluster). Otherwise, the time series with the largest 
#' contribution to the test statistics are temporarily removed (the number of time 
#' series to remove depends on the \code{rate} of removal) and \code{\link{sync.test}} 
#' is applied again. The contribution to the test statistic is assessed by the
#' WAVK test statistic calculated for each time series.
#' 
#' 
#' @param formula an object of class "\code{\link[stats]{formula}}", 
#' specifying the type of common trend 
#' for clustering the time series in a \eqn{T} by \eqn{N} matrix of time series 
#' (time series in columns) which is passed to \code{\link{sync.test}}. 
#' Variable \eqn{t} should be used to specify the form 
#' of the trend, where \eqn{t} is specified within the function automatically as a 
#' regular sequence of length \eqn{T} on the interval (0,1]. See \code{Examples}.
#' @param rate rate of removal of time series. Default is 1 (i.e., if hypothesis 
#' of synchronism is rejected one time series is removed at a time to re-test the 
#' remaining time series). Integer values above 1 are treated as number of time 
#' series to be removed. Values from 0 to 1 are treated as a fraction of 
#' time series to be removed.
#' @param alpha significance level for testing hypothesis of a common trend 
#' (using \code{\link{sync.test}}) of the parametric form specified in \code{formula}.
#' @param ... arguments to be passed to \code{\link{sync.test}}, for example, 
#' number of bootstrap replications (\code{B}).
#' 
#' 
#' @return A list with the elements:
#' \item{cluster}{an integer vector indicating the cluster to which each time series is 
#' allocated. A label \code{'0'} is assigned to time series which do not have a common trend 
#' with other time series (that is, all time series labeled with \code{'0'} are separate 
#' one-element clusters).}
#' \item{elements}{a list with names of the time series in each cluster.}
#' 
#' The further elements combine results of \code{\link{sync.test}} for each cluster with
#' at least two elements (that is, single-element clusters labeled with 
#' \code{'0'} are excluded):
#' \item{estimate}{a list with common parametric trend estimates obtained by 
#' \code{\link{sync.test}} for each cluster. 
#' The length of this list is \code{max(cluster)}.}
#' \item{pval}{a list of \eqn{p}-values of \code{\link{sync.test}} for each cluster.
#' The length of this list is \code{max(cluster)}.}
#' \item{statistic}{a list with values of \code{\link{sync.test}} test statistic for 
#' each cluster. The length of this list is \code{max(cluster)}.}
#' \item{ar_order}{a list of AR filter orders used in \code{\link{sync.test}} for each 
#' time series. The results are grouped by cluster in the list of length \code{max(cluster)}.}
#' \item{window_used}{a list of local windows used in \code{\link{sync.test}} for each 
#' time series. The results are grouped by cluster in the list of length \code{max(cluster)}.}
#' \item{all_considered_windows}{a list of all windows considered in 
#' \code{\link{sync.test}} and corresponding test results, for each cluster. 
#' The length of this list is \code{max(cluster)}.}
#' \item{WAVK_obs}{a list of WAVK test statistics obtained in \code{\link{sync.test}} 
#' for each time series. 
#' The results are grouped by cluster in the list of length \code{max(cluster)}.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{BICC}}, \code{\link{DR}}, \code{\link{sync.test}}
#' 
#' @keywords cluster trend synchrony
#' 
#' @author Srishti Vishwakarma, Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' \dontrun{
#' ## Simulate 4 autoregressive time series, 
#' ## 3 having a linear trend and 1 without a trend:
#' set.seed(123)
#' T = 100 #length of time series
#' N = 4 #number of time series
#' X = sapply(1:N, function(x) arima.sim(n = T, 
#'            list(order = c(1, 0, 0), ar = c(0.6))))
#' X[,1] <- 5 * (1:T)/T + X[,1]
#' plot.ts(X)
#' 
#' # Finding clusters with common linear trends:
#' LinTrend <- sync.cluster(X ~ t) 
#'   
#' ## Sample Output:
#' ##[1] "Cluster labels:"
#' ##[1] 0 1 1 1
#' ##[1] "Number of single-element clusters (labeled with '0'): 1"
#' 
#' ## plotting the time series of the cluster obtained
#' for(i in 1:max(LinTrend$cluster)) {
#'     plot.ts(X[, LinTrend$cluster == i], 
#'             main = paste("Cluster", i))
#' }
#' 
#' 
#' ## Simulating 7 autoregressive time series, 
#' ## where first 4 time series have a linear trend added 
#' set.seed(234)
#' T = 100 #length of time series
#' a <- sapply(1:4, function(x) -10 + 0.1 * (1:T) + 
#'             arima.sim(n = T, list(order = c(1, 0, 0), ar = c(0.6))))
#' b <- sapply(1:3, function(x) arima.sim(n = T, 
#'             list(order = c(1, 0, 0), ar = c(0.6))))
#' Y <- cbind(a, b)
#' plot.ts(Y)
#' 
#' ## Clustering based on linear trend with rate of removal = 2 
#' # and confidence level for the synchronism test 90%
#' LinTrend7 <- sync.cluster(Y ~ t, rate = 2, alpha = 0.1, B = 99)
#'    
#' ## Sample output:
#' ##[1] "Cluster labels:"
#' ##[1] 1 1 1 0 2 0 2
#' ##[1] "Number of single-element clusters (labeled with '0'): 2"
#' }
#' 
sync.cluster <- function(formula, rate = 1, alpha = 0.05, ...) 
{
    ## separating formula to find the time series
    frml <- deparse(substitute(formula))
    splt <- strsplit(frml, "~")[[1]]
    DNAME <- splt[1]
    sh <- splt[2]
    Y <- as.data.frame(eval(parse(text = DNAME)))
    OrigNames <- colnames(Y)
    # assigning column names
    N <- ncol(Y)
    colnames(Y) <- 1:N 
    # Storing the final cluster labels
    cluster <- rep(NA, N)
    # Storing outputs of sync.test for each cluster
    sync.common_trend_estimates <- 
        sync.p.value <- 
        sync.statistic <- 
        sync.ar.order_used <- 
        sync.Window_used <- 
        sync.all_considered_windows <- 
        sync.wavk_obs <- list()
    Ystar <- Y
    Ystar_names <- colnames(Ystar)
    # index for clusters
    K = 1
    while (sum(is.na(cluster)) > 1) { #while 2 or more time series are unclustered
        # testing sync.test() in the input matrix
        SyncResults <- do.call(sync.test, args = list(as.formula(paste("Ystar", "~", sh)), ...))
        # if we fail to reject the H0 of synchronism
        if (SyncResults$p.value >= alpha)
        {
            # finding common series
            j <- as.numeric(intersect(colnames(Y), colnames(Ystar)))
            # Assigning the cluster label
            cluster[j] <- K
            Ystar <- Y[, is.na(cluster)]
            Ystar_names <- colnames(Ystar)
            # Storing the output of sync.test()
            sync.common_trend_estimates[[K]] = SyncResults$estimate$common_trend_estimates
            sync.p.value[[K]] = SyncResults$p.value
            sync.statistic[[K]] = SyncResults$statistic
            sync.ar.order_used[[K]] = SyncResults$estimate$ar.order_used
            sync.Window_used[[K]] = SyncResults$estimate$Window_used
            sync.all_considered_windows[[K]] = SyncResults$estimate$all_considered_windows
            sync.wavk_obs[[K]] = SyncResults$estimate$wavk_obs
            K = K + 1
        } else {#the H0 of trend synchronism is rejected
            # Extracting local factor statistics
            WAVKResults <- abs(SyncResults$estimate$wavk_obs)
            # Sorting the WAVK result
            WAVKtmp <- order(WAVKResults, decreasing = TRUE)
            if (rate >= 1) {
                nRM <- round(rate)
            } else {
                nRM <- round(rate*length(WAVKResults))
            }
            # Make sure at least 1 time series should be removed, but less than 
            # the number of time series remaining in Ystar
            if (nRM >= ncol(Ystar) || nRM == 0) {
                nRM <- 1 
            }
            # Removing the time series as per the rate
            Ystar <- Ystar[, -WAVKtmp[1:nRM]]
            Ystar_names <- Ystar_names[-WAVKtmp[1:nRM]]
        }
        if (is.vector(Ystar)) {
            j <- as.numeric(intersect(colnames(Y), Ystar_names))
            # assigning the cluster number to cluster label variable
            cluster[j] <- 0
            Ystar <- Y[, is.na(cluster)]
            Ystar_names <- colnames(Ystar)
        }
    } #end of the 'while' loop (means at least N-1 time series are clustered)
    # If there is one time series left unclustered, assing in to a separate cluster
    cluster[is.na(cluster)] <- 0
    # Clusters themselves
    elements <- tapply(OrigNames, cluster, c)
    names(elements)[names(elements) == "0"] <- "Time series that each formed a separate cluster"
    ## Final clustering results:
    print("Cluster labels:")
    print(cluster)
    print(paste("Number of single-element clusters (labeled with '0'):", 
                sum(cluster == 0)))
    return(invisible(structure(list(cluster = cluster, 
                                    elements = elements, 
                                    estimate = sync.common_trend_estimates, 
                                    pval = sync.p.value, 
                                    statistic = sync.statistic,  
                                    ar_order = sync.ar.order_used,
                                    window_used = sync.Window_used,
                                    all_considered_windows = sync.all_considered_windows,
                                    WAVK_obs = sync.wavk_obs))))
}
