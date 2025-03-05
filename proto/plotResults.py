import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

vector = np.loadtxt("results.txt")
dim = int(sqrt(vector.size))
quad1 = vector.reshape(dim, dim)
quad2 = np.fliplr(quad1)
quad3 = np.flipud(quad1)
quad4 = np.fliplr(quad3)
matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))
plt.imshow(matrix, cmap='gray_r', interpolation='nearest')
plt.colorbar()
plt.show()