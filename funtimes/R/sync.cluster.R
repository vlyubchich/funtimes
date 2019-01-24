sync.cluster <- function(formula, rate = 1, alpha = 0.05, ...) {
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
    
    # initializing variables according to the algorithm
    Y_star <<- Y
    
    # number of columns in a matrix
    N <- ncol(Y_star)
    # number of rows in a matrix
    nrows <- nrow(Y_star)
    # index for clusters
    K = 1
    # cluster labels
    L = rep(NaN,N)
    
    # assigning column names
    colnames(Y_star) <- 1:N
    colnames(Y) <- 1:N
    
    while (!is.null(ncol(Y))) {
        if (ncol(Y_star) == 0 || is.null(ncol(Y_star))) {break}
        
        # synchronism test on Ystar
        SyncResults <- do.call(sync.test,args=list(as.formula(paste("Y_star","~",sh)), ...))
        
        # if we fail to reject the Null Hypothesis
        if (SyncResults$p.value>= alpha)
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
            
            K = K+1
            
        } else {
            # extracting local factor statistics
            WAVKResults <- abs(SyncResults$estimate$wavk_obs)
            #sorting the WAVK result
            WAVKtmp <- sort(WAVKResults)
            # rate of removal of time series
            nRM <- max(1, round(rate*length(WAVKResults)))
            WAVKtmp.rmv <- WAVKtmp[(length(WAVKtmp)-nRM+1):length(WAVKtmp)]
            # removing the time series as per the rate
            Y_star <- Y_star[, !(WAVKResults == WAVKtmp.rmv)]
        }
        if (is.vector(Y_star)) {
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
    Lfinal <- L
    clus_col.Idx <- sapply(1:max(unique(L[!is.nan(L)])), function(x) which(L == x)) 
    clus_col.NoBind <- which(L == 0)
    print(list(Clusters=table(Lfinal)))
    invisible(structure(list(Clusters = Lfinal, Column.Index = clus_col.Idx, Pval = sync.pval.Lst, 
                             TestStatistics = sync.Teststat.Lst, Estimate = sync.stat.Est.Lst, AROrder = sync.ar_order.Lst,
                             WindowUsed = sync.window_used.Lst,
                             allConsideredWindow = sync.all_consideredWindow.Lst, WAVKobs = sync.wavk_obs.Lst)))
    # removing Y_star 
    rm(Y_star, envir = parent.frame())
}
