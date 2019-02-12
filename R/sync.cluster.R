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
#' (i.e., assigned to a cluster). Othewise, the time series with the largest 
#' contribution to the test statistics are temporarily removed (the number of time 
#' series to remove depends on the \code{rate} of removal) and \code{\link{sync.test}} 
#' is applied again. The contribution to the test statistic is assessed by the
#' WAVK test statistic calculated for each time series.
#' 
#' 
#' @param formula an object of class "\code{\link[stats]{formula}}", 
#' specifying the type of common trend 
#' for clustering the time series in a \eqn{T} by \eqn{N} matrix of time series 
#' (time series in columns). It is passed to \code{\link{sync.test}}. 
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
#' @return A list with the elements: #SL: need to make the output more user-friendly. Where are the clustering results (cluster assignments or so-called cluster labels)? #SL: Also, please, use similar notations to other functions.
#' \item{Clusters}{number of clusters obtained.}
#' \item{Column.Index}{index of columns of time series (based on the main matrix) in each cluster.}
#' \item{Estimate}{parametric trend estimates of clusters obtained.}
#' \item{Pval}{\code{p}-value of the \code{sync.test}.}
#' \item{Statistics}{value of \code{sync.test} test statistics.}
#' \item{ar.order}{AR filter order.}
#' \item{window_used}{window used for the \code{sync.test} test.}
#' \item{all_considered_Windows}{different windows considered in the \code{sync.test}.}
#' \item{WAVK_obs}{WAVK test statistics for each cluster obtained from \code{sync.test}.}
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
#' T = 100
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
#' ## $`Clusters`
#' ## Lfinal
#' ## 1 
#' ## 3 
#' 
#' ## simulating seven auroregressive time series 
#' ## Three have linear trend added and four have no  trend
#' n = 50
#' nc = 7
#' Y = matrix(NA, nrow = n, ncol = nc)
#' for ( i in 1:nc){
#'     if (i < 5){
#'         Y[,i] <-  arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
#'     } else {
#'         Y[,i] <- -10 + 0.5 * (1:n) + arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
#'     }
#' }
#' 
#' plot.ts(Y)
#' ## Clustering with rate of removal = 5 and window = 15
#' \dontrun{
#'     LinTrendR5W15 <- sync.cluster(Y~t, rate = 5, Window = 15) 
#' }    
#' ## Sample output:
#' # $`Clusters`
#' # Lfinal
#' # 1 
#' # 2 
#' 
#' 
#' ## simulating five autoregressive time series to test for quadratic trend
#' ## One has no trend, while rest of the series have quadratic trend
#' n = 30
#' nc = 5
#' Y = matrix(NA, nrow = n, ncol = nc)
#' p <- 0.5
#' q <- 1:n
#' 
#' for ( i in 1:nc){
#'     if (i < 2){
#'         Y[,i] <- arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
#'     } else {
#'         Y[,i] <- -10 + p*(q+10)^2 + arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
#'     }
#'     
#' }
#' plot.ts(Y)
#' # Clustering with default rate of removal
#' \dontrun{
#'     QuadTrend <- sync.cluster(Y~poly(t,2))
#' }
#' ## Sample output:
#' ## $`Clusters`
#' ## Lfinal
#' ## 0 1 
#' ## 2 2 
#'     
sync.cluster <- function(formula, rate = 1, alpha = 0.05, ...) 
{
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
    # number of columns in a matrix
    #N <- ncol(Y_star)
    # number of rows in a matrix
    nrows <- nrow(Y_star)
    # index for clusters
    K = 1
    # cluster labels
    L = rep(NaN,N)
    
    while (!is.null(ncol(Y))) {
        #if (ncol(Y_star) == 0 || is.null(ncol(Y_star))) {break} #SL seems that this is the case of 1 TS left. Why break without assigning the last cluster?
        # synchronism test on Ystar
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
            # storing the clusters and their results
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
            #SL: This will work with a big sample in the beginning of the iterations. What if user set rate=5, and you have only 3 TS left? What if user set rate=0.1 and with 2 TS left round(rate*length(WAVKResults)) is 0 (i.e., the algorithm is stuck becase cannot remove anything)? Need to put a condition that at least 1 TS shall be removed.
            # if rate is higher than the time series left in the matrix
            # if rate becomes zero and algorithm runs infinitely
            if (nRM > ncol(Y_star) || nRM == 0){
                nRM <- 1 # atleast one time series shoud be deleted
            }
            #SL: I think this will work faster than the if() clause above: nRM <- max(1, nRM)
            
            
            # removing the time series as per the rate
            #if (rate == 1) {  #SL: This IF can be avoided For example, can sort WAVK in decreasing order, then just select/remove first nRM ones
            #    Y_star <- Y_star[, !(WAVKResults == max(WAVKResults))]
            #} else {
            WAVKtmp.rmv <- WAVKtmp$ix[(length(WAVKtmp$ix)-nRM+1):length(WAVKtmp$ix)]
            Y_star <- Y_star[, -WAVKtmp.rmv]
            #}
        }
        if (is.vector(Y_star)) {
            if (is.null(ncol(Y))) {
                break
            } else {
                # finding the position of the matching series
                for ( idx.j in 1:length(colnames(Y))){if (length(which(Y[,idx.j] == Y_star) == TRUE) == nrows){j = idx.j}}
                # Extracting the correct column name from the original Y
                clm.nm <- colnames(Y)
                # finding the index of that time series because we need to update the vector L
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
    Lfinal <- L[!L==0]
    clus_col.Idx <- sapply(1:max(unique(L[!is.nan(L)])), function(x) which(L == x)) 
    clus_col.NoBind <- which(L == 0)
    print(list(Clusters=table(Lfinal)))
    return(invisible(structure(list(Clusters = Lfinal, Column.Index = clus_col.Idx, Pval = sync.pval.Lst, 
                                    Statistics = sync.Teststat.Lst, Estimate = sync.stat.Est.Lst, ar.order = sync.ar_order.Lst,
                                    window_used = sync.window_used.Lst,
                                    all_considered_windows = sync.all_consideredWindow.Lst, WAVK_obs = sync.wavk_obs.Lst))))
    # removing Y_star 
    rm(Y_star, envir = parent.frame())
}
