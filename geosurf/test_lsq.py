
import numpy as np

a = np.array([[2,3,-1]])
b = np.array([5])

x,y,z = np.linalg.lstsq(a, b)[0]
print x,y,z

print 2*x + 3*y - z
