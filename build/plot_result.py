import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import os

files = os.listdir("output/")
for file in files:
    if not file.endswith(".txt"):
        files.remove(file)
files.insert(0, files.pop(files.index("final.txt")))

while(True):
    print("Choose which file to plot.")
    for idx, file in enumerate(files):
        print(f"\t[{idx}]\t{file}")

    target = int(input())
    target_file = files[target]
    vector = np.loadtxt("output/" + target_file)
    dim = int(sqrt(vector.size))
    quad1 = vector.reshape(dim, dim)
    quad2 = np.fliplr(quad1)
    quad3 = np.flipud(quad1)
    quad4 = np.fliplr(quad3)
    matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))
    plt.imshow(matrix, cmap='gray_r', interpolation='nearest')
    plt.colorbar()
    plt.title("Metal fractions")
    plt.show()