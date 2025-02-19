import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix



def create_rectangle_array(number_of_points, number_of_rectangle): 
    
    mesh = np.zeros((number_of_points,number_of_points))
    count = 1 
    for i in range(number_of_points):
        for j in range(number_of_points):
           mesh[i,j] = count ; count+=1 
           
    print(mesh)
    rectangles = []
    index = 1
    for i in range(1,number_of_points):
        for j in range(0,number_of_points-1):
            line = []
            line.append(mesh[i][j])
            line.append(mesh[i][j+1])
            line.append(mesh[i-1][j+1])
            line.append(mesh[i-1][j])
            rectangles.append(line)
            
    return rectangles


def create_coordinates(L,number_of_points):
    coordinates=np.zeros((number_of_points*number_of_points,2))
    h = L/(number_of_points-1)
    for i in range(number_of_points):
        for j in range(number_of_points):
            coordinates[i*number_of_points+j][0] = j*h
            coordinates[i*number_of_points+j][1] = i*h
    
    return coordinates

def bandwidth(A):
    A = csr_matrix(A)
    A = A.tocoo()  # Conversion en format COO pour accéder aux indices des non-zéros
    return max(abs(A.row - A.col))
     
     
     
L = 0.01 
outlet = 0.004    
number_of_points = 10; h = L/(number_of_points-1)
number_of_rectangle = (number_of_points-1)*(number_of_points-1)

local_matrix = [[2/3,-1/6,-1/3,-1/6],[-1/6,2/3,-1/6,-1/3],[-1/3,-1/6,2/3,-1/6],[-1/6,-1/3,-1/6,2/3]]
rectangles = create_rectangle_array(number_of_points,number_of_rectangle)
coordinates = create_coordinates(L,number_of_points)

print(rectangles)
print(coordinates)

A = np.zeros((number_of_points*number_of_points,number_of_points*number_of_points))
k_values = np.random.uniform(0.5, 1.5, size=number_of_rectangle)
for e,rectangle in enumerate(rectangles):
    k_e = k_values[e]  
    for l in range(4):
        for m in range(4):
            index_i = int(rectangle[l])
            index_j = int(rectangle[m])
            value   = local_matrix[l][m]
            A[index_i-1,index_j-1] += k_e*value


np.random.seed(42) 

F = np.zeros(number_of_points * number_of_points)

for rectangle in rectangles: 
    F_local = (h**2)*(10**7) / (4)
    
    for l in range(4):
        index_i = int(rectangle[l]) - 1  
        F[index_i] += F_local 
    

print(bandwidth(A))        
plt.imshow(A, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()

# print(F)
T = np.linalg.solve(A,F)
print(T)