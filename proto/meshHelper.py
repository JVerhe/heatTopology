import numpy as np
from optHelper import adjoint

def create_rectangle_and_mesh(number_of_points): 
    
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
    
    coordinates=np.zeros((number_of_points*number_of_points,2))
    h = L/(number_of_points-1)
    for i in range(number_of_points):
        for j in range(number_of_points):
            coordinates[i*number_of_points+j][0] = j*h
            coordinates[i*number_of_points+j][1] = i*h
    
    return coordinates


def filter_boundary_points_with_index(coordinates, L, lower=0.003, upper=0.007):

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
    h = L/(number_of_points-1)
    F = np.zeros(number_of_points * number_of_points)

    for rectangle in rectangles: 
        F_local = (h**2)*(10**7) / (4)
    
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


def objective(v,rectangles,T,K0,k_min,k_max,p):
    D = 0 
    k_values = fill_in_k(v,k_max,k_min,p)
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        T_loc = np.zeros((4,1))
        for l in range(4):
            T_loc[l][0]=T[int(rectangle[l])-1]
        
        K = np.array(K0)*k_e
        D+= ((T_loc.T)@K@T_loc)[0][0]
          
    return D


