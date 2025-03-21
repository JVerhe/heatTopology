import numpy as np
from optHelper import adjoint
import matplotlib.pyplot as plt
from mms import *

def create_rectangle_and_mesh(number_of_points): 
    ############
    ## number of points : the number of mesh point in one fourth of the mesh
    #####
    
    mesh = np.zeros((number_of_points,number_of_points))
    count = 0
    for i in range(number_of_points):
        for j in range(number_of_points):
           mesh[i,j] = count ; count+=1 
           
    rectangles = []
    for i in range(1,number_of_points):
        for j in range(0,number_of_points-1):
            line = []
            line.append(mesh[i][j])
            line.append(mesh[i][j+1])
            line.append(mesh[i-1][j+1])
            line.append(mesh[i-1][j])
            rectangles.append(line)
    
    #print(mesh) ; print(rectangles)
    return rectangles


def create_coordinates(L,number_of_points):
    ############
    ## number of points : the number of mesh point in one fourth of the mesh
    ## L : the length of the chip
    #####
    coordinates=np.zeros((number_of_points*number_of_points,2))
    h = (L/2)/(number_of_points-1)
    
    for i in range(number_of_points):
        for j in range(number_of_points):
            coordinates[i*number_of_points+j][0] = j*h
            coordinates[i*number_of_points+j][1] = i*h
    
    return coordinates


def filter_boundary_points_with_index(coordinates, L, lower=0.003, upper=0.007):
    ############
    ## coordinates : the coordinates in one fourth of the chip
    ## number of points : the number of mesh point in one fourth of the mesh
    ## L : the length of the chip
    #####
    boundary_points = []
    for idx, point in enumerate(coordinates):
        x, y = point
        if ((x == 0 or x == L) and lower <= y <= upper):
            boundary_points.append([idx, x, y])  # Stocke (index, x, y)

    return boundary_points
   
  
        
def fill_in_k(v,k_max,k_min,p): 
    k=np.zeros(len(v))
    for i in range(len(v)):
        k[i] = k_min +(k_max-k_min)*(v[i]**p)
    
    return k
    

def find_K(v,rectangles,number_of_points,local_matrix,k_min,k_max,penal):
    ############
    ## v : density vector with values between 0 and 1
    ## rectangles : list of the number of nodes for each rectangle
    ## local_matrix : the local matrix on one element (always the same)
    ## coordinates : the coordinates in one fourth of the chip
    ## number of points : the number of mesh point in one fourth of the mesh
    ## L : the length of the chip
    #####
    K = np.zeros((number_of_points*number_of_points,number_of_points*number_of_points))
    k_values = fill_in_k(v,k_max,k_min,p=penal)    
    
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        for l in range(4):
            for m in range(4):
                index_i = int(rectangle[l])
                index_j = int(rectangle[m])
                value   = local_matrix[l][m]
                K[index_i,index_j] += k_e*value
                
    return K


def find_F(rectangles,number_of_points,L):
    h = (L/2)/(number_of_points-1)
    F = np.zeros(number_of_points * number_of_points)
    for rectangle in rectangles: 
        F_local = (h**2)*(10**7)/(2)
        for l in range(4):
            index_i = int(rectangle[l])  
            F[index_i] += F_local 
    return F
            

def apply_boundary(K,F,boundary_points,T_k):
    for point in boundary_points : 
        index_of_point = point[0]
        F = F - T_k*K[:,index_of_point]
        for i in range(len(K)):
            K[i,index_of_point]=0
            K[index_of_point,i]=0
        K[index_of_point,index_of_point] =1
        F[index_of_point] = T_k
    
    return K,F

def transform_matrix(T):
    """
    This function takes a square matrix T and creates a new matrix that is (2n-1, 2m-1) in size.
    - The top-left quarter contains the original matrix.
    - The top-right quarter contains the matrix with its columns reversed, except for the last column.
    - The bottom-left quarter contains the matrix with its rows reversed, except for the last row.
    - The bottom-right quarter contains the matrix with both rows and columns reversed, without duplicating the last row or column.
    """
    n, m = T.shape
    new_matrix = np.zeros((2 * n - 1, 2 * m - 1))

    # Top-left: Original matrix
    new_matrix[:n, :m] = T

    # Top-right: Columns reversed (except last column)
    new_matrix[:n, m:] = T[:, -2::-1]

    # Bottom-left: Rows reversed (except last row)
    new_matrix[n:, :m] = T[-2::-1, :]

    # Bottom-right: Rows and columns reversed (except last row and last column)
    new_matrix[n:, m:] = T[-2::-1, -2::-1]

    return new_matrix



def objective(v,rectangles,T,K0,k_min,k_max,p):
    D = 0 
    k_values = fill_in_k(v,k_max,k_min,p)
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        T_loc = np.zeros((4,1))
        for l in range(4):
            T_loc[l][0]=T[int(rectangle[l])]
        
        K = np.array(K0)*k_e
        D+= ((T_loc.T)@K@T_loc)[0][0]
          
    return D


def load_result_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        U = np.array([float(line.strip()) for line in lines])
    return U




L = 0.01
number_of_points = 250
p = 1
T_k =  293
period = 1
k_constant = 0.5
local_matrix = [[2/3,-1/6,-1/3,-1/6],[-1/6,2/3,-1/6,-1/3],[-1/3,-1/6,2/3,-1/6],[-1/6,-1/3,-1/6,2/3]]

rectangles = create_rectangle_and_mesh(number_of_points) ; v = np.ones(len(rectangles))*(0.8/(65-0.2))
coordinates = create_coordinates(L,number_of_points)
#boundary_points = filter_boundary_points_with_index(coordinates,L)
boundary_points = filter_boundary_points_with_index_mms(coordinates,L)


mms = True

if mms==False : 
    h = (L/2)/(number_of_points-1)
    q_rectangle = create_q_rectangle_middle(rectangles,coordinates,L,number_of_points,period)
    
    #q_rectangle = create_q_rectangle_middle_k_variable(rectangles,coordinates,L,number_of_points,period)
    #K = find_K_with_k_variable(v,rectangles,number_of_points,local_matrix,coordinates,0.2,65,p)
    F = find_F_mms(rectangles,number_of_points,L,q_rectangle)
    
    K = find_K(v,rectangles,number_of_points,local_matrix,0.2,65,p)
    #F = find_F(rectangles,number_of_points,L)
    
    K,F = apply_boundary(K,F,boundary_points,T_k)
    T_matrix = np.linalg.solve(K,F)
    T_matrix_r = T_matrix.reshape((number_of_points, number_of_points))
    #T_matrix_r = transform_matrix(T_matrix_r)


if mms==True:
    
    h = (L)/(number_of_points-1)
    T_matrix = load_result_from_file("../build/output/temperature_mms.txt")
    T_matrix_r = T_matrix.reshape((number_of_points, number_of_points))


plt.figure(figsize=(6, 5))
plt.imshow(T_matrix_r, cmap='magma', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label="Température (K)")
plt.title("Distribution de la température")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.show()



if mms==True:

    ##T_true  = create_T_point(coordinates,k_constant,L,number_of_points,period)
    T_true = load_result_from_file("../build/output/correct_temperature_mms.txt")
    T_true_r = T_true.reshape((number_of_points, number_of_points))
    Err  = T_matrix_r - T_true_r

    plt.figure(figsize=(6, 5))
    plt.imshow(Err, cmap='gray_r', origin='lower', extent=[0, L, 0, L])
    plt.colorbar(label="Error (K)")
    plt.title("Error with the true solution")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.show()



    # Points de quadrature de Gauss pour un carré de référence [-1,1]x[-1,1]
    xi_q = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    eta_q = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    w_q = np.ones(4)  


    error_L2 = 0.0

    def shape_functions(xi, eta):
        N1 = (1 - xi) * (1 - eta) / 4
        N2 = (1 + xi) * (1 - eta) / 4
        N3 = (1 + xi) * (1 + eta) / 4
        N4 = (1 - xi) * (1 + eta) / 4
        return np.array([N1, N2, N3, N4])


    number_of_elements = np.sqrt(len(rectangles))
    for i in range(len(rectangles)):
        rectangle = rectangles[i]
        i1 = int(rectangle[0])
        i2 = int(rectangle[1])
        i3 = int(rectangle[2])
        i4 = int(rectangle[3])
        Uh_nodes = np.array([T_matrix[i1],T_matrix[i2],T_matrix[i3],T_matrix[i4]])
        U_nodes = np.array([T_true[i1],T_true[i2],T_true[i3],T_true[i4]])
        
        
        for xi in xi_q:
            for eta in eta_q:
                N = shape_functions(xi, eta)  
                U_q = np.dot(N, U_nodes) 
                Uh_q = np.dot(N, Uh_nodes) 
                Err_2 = (Uh_q-U_q)**2
                #Err_2 = np.abs(U_q-Uh_q)
                J = (h**2) / 9
                error_L2 += np.sum(Err_2)* J


    error_L2 = np.sqrt(error_L2)
    print("Error :" +str(error_L2))
    print("Error/h^2 :" +str(error_L2/(h**2)))
    
    print(np.log10(error_L2))
    print(np.log10(h))
    
