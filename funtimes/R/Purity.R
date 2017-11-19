purity <- function(classes, clusters) 
{
  tmp <- cbind(classes, clusters)
  tmp <- tmp[order(tmp[,1]), ]
  classes <- tmp[,1]
  clusters <- tmp[,2]
  nClasses <- length(unique(classes))
  nClusters <- length(unique(clusters))
  ClassLabels <- ClusterLabels <- ClusterSize <- numeric()
  tmp <- table(classes, clusters)
  for(i in 1:min(dim(tmp))){
    ind <- which(tmp==max(tmp), arr.ind=TRUE)
    if(nrow(ind)>1){ #several equal max values
      if(any(dim(tmp)<2)) {
        ind <- ind[1,]
      }else{
        #for each class with equal max, run and select the second max and choose minimum of that
        ind <- ind[which.min(sapply(1:nrow(ind), function(j) sort(tmp[,ind[j,2]],decreasing=TRUE)[2])), ]
      }
    }
    ClassLabels[i] <- rownames(tmp)[ind[1]]
    ClusterLabels[i] <- colnames(tmp)[ind[2]]
    ClusterSize[i] <- max(tmp)
    tmp <- tmp[-which(rownames(tmp)==ClassLabels[i]), -which(colnames(tmp)==ClusterLabels[i])]
    #If converged to only one column (cluster), tmp should be a column-vector:
    if(is.null(dim(tmp))){
      if(nClusters < nClasses){ 
        tmp <- as.matrix(tmp)
        colnames(tmp) <- unique(clusters)[!is.element(unique(clusters), ClusterLabels)]
      }else{ #tmp should be a row-vector
        tmp <- t(as.matrix(tmp))
        rownames(tmp) <- unique(classes)[!is.element(unique(classes), ClassLabels)]
        if(length(tmp)==1) colnames(tmp) <- unique(clusters)[!is.element(unique(clusters), ClusterLabels)]
      }
    }
  }
  out <- data.frame(ClassLabels, ClusterLabels, ClusterSize)
  out <- out[order(out$ClassLabels),]
  list(pur=sum(ClusterSize)/length(classes), out=out)
}
