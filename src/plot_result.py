import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import os
import re


def listFiles():
    files = os.listdir("output/")

    if len(files) == 1:
        print("No files in the output folder")
        exit()

    for file in files:
        if not file.endswith(".txt"):
            files.remove(file)

    files.remove("results.txt")
    files = sorted(files, key=lambda x: int(re.findall(r"iteration\s*(.*?)(?=\.)", x)[-1]))
    files.append("results.txt")
    # files.insert(0, files.pop(files.index("results.txt")))

    print("\n\n")
    print("Choose which file to plot. (Press ctrl+C to exit)")

    print("\t[0]\tPlot most optimal solution")
    print("\t[1]\tPlot evolution of solutions")
    print("\t[2]\tExit")

    return files

# TODO: input validation
def listIntermediateSolutions(files):

    print("\nType the indices of the graphs you want to plot separated with a space, for example: '1 2 3'\n")

    for idx, file in enumerate(files):
        print(f"\t[{idx}]\t{file}")

    idxList = list()
    input_string = input()
    string_list = input_string.split()
    idxList = [int(num) for num in string_list]

    return idxList


def plotOptimalSolution(files):
    target_file = files[0]
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
    return


def plotEvolution(idxList, files):
    size = len(idxList)
    figure, axis = plt.subplots(1, size)

    for i, fileIdx in enumerate(idxList):
        target_file = files[fileIdx]
        
        vector = np.loadtxt("output/" + target_file)
        dim = int(sqrt(vector.size))
        quad1 = vector.reshape(dim, dim)
        quad2 = np.fliplr(quad1)
        quad3 = np.flipud(quad1)
        quad4 = np.fliplr(quad3)
        matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))

        axis[i].imshow(matrix, cmap='gray_r', interpolation='nearest')

        if target_file != "results.txt":
            iterationNr = re.findall(r"iteration\s*(.*?)(?=\.)", target_file)
            axis[i].set_title("Iteration " + iterationNr[0])

        else:
            axis[i].set_title("Final result")

    plt.show()
    return


def readUserInput(files):
    try:
        target = int(input())
    except ValueError:
        print("Not a valid input.")
        return
    
    match target:
        case 0:
            plotOptimalSolution(files)
            return
        case 1:
            plotIndices = listIntermediateSolutions(files)
            plotEvolution(plotIndices, files)
            return
        case 2:
            exit()
        case _:
            print("Not a valid input.")
            return

if __name__ == "__main__":
    files = listFiles()
    while(True):
        readUserInput(files)
        files = listFiles()

