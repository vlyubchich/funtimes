DR <- function(X, method, minPts = 3, theta = 0.9, B = 500, lb = -30, ub = 10)
{
  control_para <- sapply(seq(lb,ub,by=0.5),function(x) 1.1^x)
  Nnodes <- floor(ncol(X)/2)
  
  # define local minimum selection function
  para_opt <- function(M){
    acd <- M[,2]
    max <- max(acd)
    max_idx <- which(acd == max)
    max_idx <- max_idx[length(max_idx)]
    if (max == 0) return("All ACDs are 0. This may be caused by the range of controlling parameters is either too small or too large.")
    else if(max_idx > (length(acd)-6)) return("Max ACD is found near the upper bound of controlling parameters, there is no room to check the local minimum. You might need increase the upper bound of controlling parameters.")
    else{
      for (i in (max_idx+3):(length(acd)-3)){
        if(acd[i-1] >= acd[i] & acd[i] < acd[i+1]){
          if(acd[i-2] >= acd[i] & acd[i] < acd[i+2]){
            if(acd[i-3] >= acd[i] & acd[i] < acd[i+3]) return(M[i,1])
          }
        }
      }
      return("Not finding a local minimum.")
    }
  }
  
  cat(sprintf("Total # of controlling parameter: %d \n",length(control_para)))
  cat(sprintf("The # of controlling parameters has been processed: \n"))
  
  if(method == "DBSCAN"){
    if (missing(minPts)) stop("missing parameter minPts for DBSCAN")
    else{
      
      ACD <- lapply(1:length(control_para),function(x){
        Buffer <- c()
        for(i in 1:B){
          Sindex<-sample(1:ncol(X),Nnodes,replace=F)
          sub_data_1 <- X[,Sindex]
          cl_1 <- max(dbscan::dbscan(t(sub_data_1),eps=control_para[x], minPts=minPts)$cluster)
          
          sub_data_2 <- X[,-Sindex]
          cl_2 <- max(dbscan::dbscan(t(sub_data_2),eps=control_para[x], minPts=minPts)$cluster)
          Buffer <- c(Buffer,abs(cl_1-cl_2))
        }
        print(x)
        return(list(control_para[x],mean(Buffer)))
      })
    }
    
  }else if(method == "TRUST"){
    
    ACD <- lapply(1:length(control_para),function(x){
      Buffer <- c()
      for(i in 1:B){
        
        Sindex<-sample(1:ncol(X),Nnodes,replace=F)
        sub_data_1 <- X[,Sindex]
        cl_1 <- max(CSlideCluster(sub_data_1,Delta=control_para[x],Theta=theta))
        
        sub_data_2 <- X[,-Sindex]
        cl_2 <- max(CSlideCluster(sub_data_2,Delta=control_para[x],Theta=theta))
        
        Buffer <- c(Buffer,abs(cl_1-cl_2))
      }
      print(x)
      return(list(control_para[x],mean(Buffer)))
    })
  }else stop ("Missing method")
  
  ACD_matrix <- apply(do.call('rbind',ACD),2,as.numeric)
  colnames(ACD_matrix) <- c("Controlling Parameter","ACD")
  p_opt <- para_opt(ACD_matrix)
  out <- list(p_opt, ACD_matrix)
  names(out) <- c("P_opt","ACD_matrix")
  return(out)
}


