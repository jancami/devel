import numpy as np

def powered_sum(n,k):
	result = np.arange(n+1)
	print(result)
	for loop in range(n):
		sum=0
		for subloop in range(1,k): 
			sum = sum+subloop**float(loop)
			print(sum)
		result[loop]=sum
		print(loop, result)
	print(result)

if __name__ == "__main__":
	a=powered_sum(4,3)