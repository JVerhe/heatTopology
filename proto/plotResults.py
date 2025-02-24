import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

vector = np.loadtxt("results.txt")
dim = int(sqrt(vector.size))
matrix = vector.reshape(dim, dim)
plt.imshow(matrix, cmap='gray_r', interpolation='nearest')
plt.colorbar()
plt.show()