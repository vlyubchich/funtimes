import sys
sys.setrecursionlimit(100000)
import numpy as np

def _hist_find_new(l, MAX_STEP, n_bin):
    
    srch,perctl = np.histogram(l, bins=n_bin)
    
    key_perctl = None
    
    for i in np.arange(len(srch)-1-MAX_STEP,MAX_STEP+1,-1):
        
        tmp_bool = []
        for step in range(1, MAX_STEP+1):
            #if srch[i] < srch[i+step] and srch[i-1] < srch[i+step] and max(srch[i], srch[i-1]) <= srch[i+step]*(1-stepness): 
            if srch[i] < srch[i-step] and srch[i] < srch[i+step]:# and (srch[i] <= srch[i-step-1] or srch[i] <= srch[i-step-2]):
                tmp_bool.append(True)
            else:
                tmp_bool.append(False)
                
        tmp_bool = np.array(tmp_bool)
        num_ = np.where(tmp_bool == True)[0].tolist()

        if len(num_) == MAX_STEP:

            slope = float(srch[i+MAX_STEP] - srch[i])#/float(perctl[1] - perctl[0])
            key_perctl = perctl[i]
            break
    
    
    idx = np.where(l > key_perctl)[0].tolist()
    return idx



def cal_adjM_cutOff(xxDist, MAX_STEP, num_bin):

    adj = []
    for k in range(xxDist.shape[0]):

        temp = [k]
        
        l = np.array(xxDist[k,])
        
        idx = _hist_find_new(l, MAX_STEP, num_bin)  
        #idx = np.where(l > key_perctl)[0].tolist()
        
        if k in idx:
            idx.remove(k)

        temp.extend(idx)

        adj.append(temp)
        
    return adj

#clustering algorithm
def dfs(n, cl, adj, label):

	if cl[n] == -1:

		cl[n] = label

		if len(adj[n])-1 !=0:
			
			for j in range(1, len(adj[n])):

				dfs(adj[n][j], cl, adj, label)

	return

def clustering_(adj):
    n = len(adj)
    cl = [-1]*n
    label = 0

    for i in range(n):

        if cl[i] == -1:


            if len(adj[i]) == 1:
                cl[i] = 0

            else:

                label = label+1
                dfs(i, cl, adj, label)
                
    return cl


def CRAD(data, MAX_STEP, num_bin, cov):

    data = np.array(data)
    cov = np.array(cov)
    try:
        VI = np.linalg.inv(cov)
    except:
        VI = np.linalg.pinv(cov)
    
    
    n = data.shape[0]
    xxDist = np.zeros((n,n))-1

    
    for i in range(n):
        x1 = np.array(data[i,:],ndmin=2)
        for j in range(i,n):
            
            x2 = np.array(data[j,:],ndmin=2)
            c = np.array(x2-x1,ndmin=2)
            dist = 1.0/ (1+np.dot(np.dot(c,VI),c.T))
            
            xxDist[i][j] = dist
            xxDist[j][i] = dist

    adj = cal_adjM_cutOff(xxDist, MAX_STEP, num_bin)

    cl = clustering_(adj)

    return cl
