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
#' regular sequence of length \eqn{T} on the interval (0,1]. See `Examples'.
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
#'  #SL: need to make the output more user-friendly. Where are the clustering results (cluster assignments or so-called cluster labels)? #SL: Also, please, use similar notations to other functions.
#' \item{clusters}{total number of clusters obtained.}
#' \item{cluster_label}{cluster labels of each time series (each label is the cluster number of 
#' time series). A label zero represents the time series which does not have a common trend with
#' other time series. Labels other than zero are the clusters obtained.}
#' \item{column_index}{stores index of columns of time series in the main matrix corresponding to each cluster. This can be used for
#' plotting the time series of clustered time series, for example.}
#' \item{estimate}{parametric trend estimates obtained for each cluster.}
#' \item{pval}{\code{p}-value of \code{sync.test} for each cluster.}
#' \item{statistics}{value of \code{sync.test} test statistics for each cluster.}
#' \item{ar_order}{AR filter order used in \code{sync.test} for each cluster.}
#' \item{window_used}{window used in the \code{sync.test} test for each cluster.}
#' \item{all_considered_Windows}{different windows considered in \code{sync.test} for each cluster.}
#' \item{WAVK_obs}{WAVK test statistics for each cluster obtained from \code{sync.test} for each cluster.}
#' In case of more than one cluster, each output will be in 
#' the same sequence as of unique cluster labels except zero.
#' 
#' 
#' @references
#' \insertAllCited{}
#' 
#' @keywords cluster trend synchrony
#' 
#' @author Srishti Vishwakarma, Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' ## Simulate four autoregressive time series, 
#' ## three having a linear trend and one without a trend:
#' set.seed(123)
#' T = 30
#' N = 4
#' X = sapply(1:N, function(x) arima.sim(n = T + 100, list(order = c(1, 0, 0), ar = c(0.6)))[-c(1:100)])
#' X[,1] <- 5 * (1:T)/T + X[,1]
#' plot.ts(X)
#' 
#' # Finding clusters with default parameters
#' \dontrun{
#'     LinTrend <- sync.cluster(X ~ t) 
#' }    
#' ## Sample Output:
#' ## "Number of Clusters obtained: 1"
#' ## "Cluster Labels"
#' ## 1 1 0 0
#' 
#' ## plotting the time series of the cluster obtained
#' idx = LinTrend$column_index
#' plot.ts(X[,idx])
#' 
#' ## simulating seven autoregressive time series 
#' ## Four time series have linear trend added 
#' set.seed(123)
#' T = 30
#' N = 10
#' X = matrix(NA, nrow = T, ncol = N)
#' X[,1:5] <-  sapply(1:5, function(x) arima.sim(n = T + 100, list(order = c(1, 0, 0), ar = c(0.6)))[-c(1:100)])
#' X[,6:N] <- sapply(6:N, function(x) -10 + 0.5 * (1:T) + arima.sim(n = T + 100, list(order = c(1, 0, 0), ar = c(0.6)))[-c(1:100)])

#' plot.ts(X)
#' 
#' ## Clustering with rate of removal = 5, and window = 15
#' \dontrun{
#'     LinTrendR5W15 <- sync.cluster(X~t, rate = 2, Window = 15) 
#' }    
#' ## Sample output:
#' ## "Number of Clusters obtained: 2"
#' ## "Cluster Labels"
#' ##   2 1 1 2 1 1 1 1 1 1
#' 
#' ## simulating five autoregressive time series to test for quadratic trend
#' ## One has no trend, while rest of the series have quadratic trend
#' set.seed(123)
#' T = 30
#' N = 5
#' X = matrix(NA, nrow = T, ncol = N)
#' p <- 0.5
#' q <- 1:T
#' 
#' X[,1:2] <-  sapply(1:2, function(x) arima.sim(n = T + 100, list(order = c(1, 0, 0), ar = c(0.6)))[-c(1:100)])
#' X[,3:N] <- sapply(3:N, function(x) -10 + p*(q+10)^2 +arima.sim(n = T + 100, list(order = c(1, 0, 0), ar = c(0.6)))[-c(1:100)])
#'   
#' plot.ts(X)
#' # Clustering with default rate of removal
#' \dontrun{
#'     QuadTrend <- sync.cluster(X~poly(t,2))
#' }
#' ## Sample output:
#' ## "Number of Clusters obtained: 1"
#' ## "Cluster Labels"
#' ## 1 1 0 1 1  
#' 
#' 
sync.cluster <- function(formula, rate = 1, alpha = 0.05, ...) 
{
    require(funtimes)
    # Storing the final list of clusters 
    Lfinal <- list() 
    clus_col.Idx <- list() # Storing the index of columns in cluster
    # Storing values for each category
    sync.pval.Lst = list()
    sync.stat.Est.Lst = list()
    sync.Teststat.Lst = list()
    sync.ar_order.Lst = list()
    sync.window_used.Lst = list()
    sync.all_consideredWindow.Lst = list()
    sync.wavk_obs.Lst = list()
    ## separating formula to find the time series
    frml <- deparse(substitute(formula))
    splt <- strsplit(frml,"~")[[1]]
    DNAME <- splt[1]
    sh <- splt[2]
    Y <- as.data.frame(eval(parse(text = DNAME))) # Reading data
    
    # assigning column names
    N <- ncol(Y)
    colnames(Y) <- 1:N
    
    # initializing variables according to the algorithm
    Y_star <<- Y
    
    # number of rows in a matix
    nrows <- nrow(Y_star)
    # index for clusters
    K = 1
    # cluster labels
    L = rep(NaN,N)
    
    while (!is.null(ncol(Y))) {
        
        ## removing time series with certain rate may lead to empty matrix for few cases
        ## checking it first
        if (length(Y_star) == 0 ) {break}
        
        # testing sync.test() in the input matrix
        SyncResults <- do.call(sync.test, args = list(as.formula(paste("Y_star", "~", sh)), ...))
        # if we fail to reject the Null Hypothesis
        if (SyncResults$p.value >= alpha)
        {
            # finding common series
            j = intersect(colnames(Y),colnames(Y_star))
            # assigning the cluster number to cluster label variable
            j1 = as.numeric(j)
            L[j1] = K
            #removing the series from the cluster
            Y <- Y[ , !(names(Y) %in% j)]
            # updating Y star
            Y_star <- Y
            N <- ncol(Y_star)
            # storing the estimates of parametric trend in clustered time series 
            # along with other test statistics from sync.test()
            sync.stat.Est.Lst[[K]] = SyncResults$estimate$common_trend_estimates
            sync.pval.Lst[[K]] = SyncResults$p.value
            sync.Teststat.Lst[[K]] = SyncResults$statistic
            sync.ar_order.Lst[[K]] = SyncResults$estimate$ar.order_used
            sync.window_used.Lst[[K]] = SyncResults$estimate$Window_used
            sync.all_consideredWindow.Lst[[K]] = SyncResults$estimate$all_considered_windows
            sync.wavk_obs.Lst[[K]] = SyncResults$estimate$wavk_obs
            K = K + 1
        } else {
            # extracting local factor statistics
            WAVKResults <- abs(SyncResults$estimate$wavk_obs)
            #sorting the WAVK result
            WAVKtmp <- sort(WAVKResults, index.return = TRUE)
            if (rate >= 1) {
                nRM <- round(rate)
            } else {
                nRM <- round(rate*length(WAVKResults))
            } 
           
            # if the rate of removal becomes greater than the number of columns or becomes zero in the testing data
            # then rate becomes equal to 1.
            # atleast 1 time series should be removed 
            if (nRM > ncol(Y_star) || nRM == 0){
                nRM <- 1 
            }
            #SL: This will work with a big sample in the beginning of the iterations.
            #What if user set rate=5, and you have only 3 TS left?
            #What if user set rate=0.1 and with 2 TS left round(rate*length(WAVKResults)) 
            #is 0 (i.e., the algorithm is stuck becase cannot remove anything)? 
            #Need to put a condition that at least 1 TS shall be removed.
            # if rate is higher than the time series left in the matrix
            # if rate becomes zero and algorithm runs infinitely
       
            #SL: I think this will work faster than the if() clause above: nRM <- max(1, nRM)
            
            
            # removing the time series as per the rate
            WAVKtmp.rmv <- WAVKtmp$ix[(length(WAVKtmp$ix)-nRM+1):length(WAVKtmp$ix)]
            Y_star <- Y_star[, -WAVKtmp.rmv]
    
        }
        if (is.vector(Y_star)) {
            if (is.null(ncol(Y))) {
                break
            } else {
                # finding the position of the matching series
                for ( idx.j in 1:length(colnames(Y))){if (length(which(Y[,idx.j] == Y_star) == TRUE) == nrows){j = idx.j}}
                # Extracting the correct column name from the original Y
                clm.nm <- colnames(Y)
                # finding the index of that time series because we need to update the vector L (cluster label)
                j1 = as.numeric(clm.nm[j])
                # updating the L with correct index of time series
                L[j1] = 0 # zero denotes that time series is not joined by any other series
                # removing the matching series
                Y <- Y[ , -j]
                Y_star <- Y
                N <- ncol(Y_star)
                
            }
        }
      
        Y_star <<- Y_star
        
    }
    Lfinal <- L 
    ## The very last time series in the loop will be left out without being alloted to cluster
    Lfinal[is.nan(Lfinal)] = 0
    clus_col.Idx <- sapply(1:max(unique(Lfinal)), function(x) which(Lfinal == x)) 
    ## final cluster variable
    print(paste("Number of Clusters obtained: ", max(Lfinal), sep = ""))
    ## Cluster assginment
    print("Cluster Labels")
    print(Lfinal)
    
    return(invisible(structure(list(Clusters = max(Lfinal), cluster_label = Lfinal, column_index = clus_col.Idx, 
                                    estimate = sync.stat.Est.Lst, pval = sync.pval.Lst, 
                                    statistics = sync.Teststat.Lst,  ar_order = sync.ar_order.Lst,
                                    window_used = sync.window_used.Lst,
                                    all_considered_windows = sync.all_consideredWindow.Lst, WAVK_obs = sync.wavk_obs.Lst))))
    # removing Y_star 
    rm(Y_star, envir = parent.frame())
}
