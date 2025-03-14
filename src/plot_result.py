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

    nonSortedFiles = ["objective_values.txt", "temperature.txt", "density.txt"]
    
    for f in nonSortedFiles:
        files.remove(f)
    files = sorted(files, key=lambda x: int(re.findall(r"iteration\s*(.*?)(?=\.)", x)[-1]))
    files.extend(nonSortedFiles)

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
    target_file = files[-1]
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


def plotCollection(idxList, files):
    size = len(idxList)
    _, axis = plt.subplots(1, size, squeeze=False)

    for i, fileIdx in enumerate(idxList):
        target_file = files[fileIdx]
        
        nonSortedFiles = ["objective_values.txt", "temperature.txt", "density.txt"]
        if target_file not in nonSortedFiles:

            vector = np.loadtxt("output/" + target_file)
            dim = int(sqrt(vector.size))
            quad1 = vector.reshape(dim, dim)
            quad2 = np.fliplr(quad1)
            quad3 = np.flipud(quad1)
            quad4 = np.fliplr(quad3)
            matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))

            axis[0,i].imshow(matrix, cmap='gray_r', interpolation='nearest')
            iterationNr = re.findall(r"iteration\s*(.*?)(?=\.)", target_file)
            axis[0,i].set_title("Iteration " + iterationNr[0])

        elif target_file == "density.txt":
            vector = np.loadtxt("output/" + target_file)
            dim = int(sqrt(vector.size))
            quad1 = vector.reshape(dim, dim)
            quad2 = np.fliplr(quad1)
            quad3 = np.flipud(quad1)
            quad4 = np.fliplr(quad3)
            matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))

            axis[0,i].imshow(matrix, cmap='gray_r', interpolation='nearest')
            axis[0,i].set_title("Final Metal fraction")

        elif target_file == "objective_values.txt":
            vector = np.loadtxt("output/" + target_file)
            axis[0, i].plot(vector)
            axis[0, i].set_yscale('log')
            axis[0, i].set_title("Cost function")
            axis[0, i].set_xlabel("Iteration")

        # TODO
        elif target_file == "temperature.txt":
            vector = np.loadtxt("output/" + target_file)
            axis[0, i].plot(vector)
            axis[0, i].set_xlabel("Iteration")
            axis[0,i].set_title("Maximal Temperature (K)")

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
            plotCollection(plotIndices, files)
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

