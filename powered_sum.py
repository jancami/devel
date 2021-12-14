import numpy as np
import time

def powered_sum4(n,k):
    array1 = np.tile(np.arange(float(1),k+1),(n+1,1))
    powers = np.reshape(np.repeat(np.arange(n+1), k), (n+1, k))
    result = np.power(array1,powers)
    result = result.sum(axis=1)
    return(result)

def powered_sum3(n,k):
    array1 = np.tile(np.arange(float(1),k+1),(n+1,1))
    powers = np.reshape(np.repeat(np.arange(n+1), k), (n+1, k))
    result = (array1**powers).sum(axis=1)
    return(result)

def powered_sum2(n,k):
    result = np.zeros(n+1)
    #print(result)
    for loop in range(n+1):
        result[loop] = ((np.arange(1,k+1))**float(loop)).sum()
    return(result)

def powered_sum(n,k):
    result = np.zeros(n+1)
    #print(result)
    for loop in range(n+1):
        sum=0.
        for subloop in range(1,k+1): 
            sum = sum+float(subloop)**float(loop)
            #print(loop,subloop,sum)
        result[loop]=sum
        #print(loop, result)
    #print(result)
    return(result)

if __name__ == "__main__":
    n=3
    k=3
    t0=time.time()
    a1=powered_sum(n,k)
    t1 = time.time()
    a2=powered_sum2(n,k)
    t2 = time.time()
    a3=powered_sum3(n,k)
    t3 = time.time()
    a4=powered_sum4(n,k)
    t4 = time.time()
    
    print(a1, a2, a3, a4)
    print(t1-t0, t2-t1, t3-t2, t4-t3)