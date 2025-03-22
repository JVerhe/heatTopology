import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns


def listFiles():
    files = os.listdir("output/")
    
    if len(files) == 1:
        print("No files in the output folder")
        exit()
    
    files = [file for file in files if file.endswith(".txt")]
    
    nonSortedFiles = ["objective_values.txt", "temperature.txt", "density.txt"]
    sortedFiles = [f for f in files if f not in nonSortedFiles]
    
    def extract_iteration(filename):
        match = re.findall(r"iteration\s*(\d+)", filename)
        return int(match[-1]) if match else float('inf')
    
    sortedFiles = sorted(sortedFiles, key=extract_iteration)
    files = sortedFiles + nonSortedFiles
    
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


def plotEvolution(idxList, files):
    size = len(idxList)
    _, axis = plt.subplots(1, size, squeeze=False)

    for i, fileIdx in enumerate(idxList):
        target_file = files[fileIdx]
        
        nonSortedFiles = ["objective_values.txt", "temperature.txt", "density.txt"]
        if target_file not in nonSortedFiles:

            vector = np.loadtxt("output/" + target_file)
            dim = int(sqrt(vector.size))
            # quad1 = vector.reshape(dim, dim)
            # quad2 = np.fliplr(quad1)
            # quad3 = np.flipud(quad1)
            # quad4 = np.fliplr(quad3)
            # matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))
            
            matrix = vector.reshape(dim,dim)

            extent = [0, 0.01, 0, 0.01] 
            img = axis[0,i].imshow(matrix, cmap='magma', interpolation='nearest',extent=extent)
            iterationNr = re.findall(r"iteration\s*(.*?)(?=\.)", target_file)
            plt.colorbar(img, ax=axis[0, i])


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
            sns.set_style("whitegrid")
            
            vector = np.loadtxt("output/objective_values.txt")
            axis[0, i].plot(vector, color='royalblue', linestyle='-', marker='o', markersize=4, alpha=0.7)
            axis[0, i].set_yscale('log')

            axis[0, i].set_title("Cost Function (Log Scale)", fontsize=14, fontweight='bold')
            axis[0, i].set_xlabel("Iteration", fontsize=12)
            axis[0, i].set_ylabel("Objective Value", fontsize=12)

            axis[0, i].grid(which="both", linestyle="--", linewidth=0.5, alpha=0.7)
            
            plt.style.use('default')
            

        # TODO
        elif target_file == "temperature.txt":
            vector = np.loadtxt("output/" + target_file)
            dim = int(sqrt(vector.size))
            quad1 = vector.reshape(dim, dim)
            quad2 = np.fliplr(quad1)
            quad3 = np.flipud(quad1)
            quad4 = np.fliplr(quad3)
            matrix = np.vstack((np.hstack((quad1, quad2)), np.hstack((quad3, quad4))))
            
            img=axis[0,i].imshow(matrix, cmap='magma', interpolation='nearest')
            axis[0,i].set_title("Temperature (K)")
            plt.colorbar(img, ax=axis[0, i])

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

