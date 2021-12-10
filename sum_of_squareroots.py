import numpy as np

def sum_of_squareroots(x,y):
	x_ok = type(x) == int or float
	y_ok = type(y) == int or float
	input_ok = x_ok and y_ok
	if input_ok: 
		print("Input OK -- continuing....")
		x_pos = x>0
		y_pos = y>0
		if x_pos and y_pos: retval = np.sqrt(x)+np.sqrt(y)
		if (x*y)<0: retval = abs(np.sqrt(abs(x)) - np.sqrt(abs(y)))
		if x_pos == False and y_pos == False: retval=1.0/(np.sqrt(abs(x))+np.sqrt(abs(y)))
		return retval


if __name__ == "__main__":
	print("Main")
	print(sum_of_squareroots(16,25))
	print(sum_of_squareroots(16,-25))
	print(sum_of_squareroots(-16,-25))