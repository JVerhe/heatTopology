import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.optimize import minimize



L = 0.01 
outlet = 0.004    
number_of_points =40; h = L/(number_of_points-1)
number_of_rectangle = (number_of_points-1)*(number_of_points-1)
p = 2
k_min = 0.2
k_max = 65
T_k =  293

local_matrix = [[2/3,-1/6,-1/3,-1/6],[-1/6,2/3,-1/6,-1/3],[-1/3,-1/6,2/3,-1/6],[-1/6,-1/3,-1/6,2/3]]


def create_rectangle_and_mesh(number_of_points): 
    
    mesh = np.zeros((number_of_points,number_of_points))
    count = 1 
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
            
    return rectangles



def create_coordinates(L,number_of_points):
    
    coordinates=np.zeros((number_of_points*number_of_points,2))
    h = L/(number_of_points-1)
    for i in range(number_of_points):
        for j in range(number_of_points):
            coordinates[i*number_of_points+j][0] = j*h
            coordinates[i*number_of_points+j][1] = i*h
    
    return coordinates


def filter_boundary_points_with_index(coordinates, L, number_of_points, lower=0.003, upper=0.007):

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
        

rectangles = create_rectangle_and_mesh(number_of_points)
coordinates = create_coordinates(L,number_of_points)
boundary_points = filter_boundary_points_with_index(coordinates,0.01,number_of_points)

    
def find_K(v,number_of_points):
    
    K = np.zeros((number_of_points*number_of_points,number_of_points*number_of_points))
    fill_in_k(v,k_max,k_min,2)
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        for l in range(4):
            for m in range(4):
                index_i = int(rectangle[l])
                index_j = int(rectangle[m])
                value   = local_matrix[l][m]
                K[index_i-1,index_j-1] += k_e*value
                
    return K


def find_F():
    
    F = np.zeros(number_of_points * number_of_points)

    for rectangle in rectangles: 
        F_local = (h**2)*(10**7) / (4)
    
        for l in range(4):
            index_i = int(rectangle[l]) - 1  
            F[index_i] += F_local 
    
    return F
            

def apply_boundary(K,F):
    
    for point in boundary_points : 
        index_of_point = point[0]
        F = F - T_k*K[:,index_of_point]
        for i in range(len(K)):
            K[i,index_of_point]=0
            K[index_of_point,i]=0
        K[index_of_point,index_of_point] =1
        F[index_of_point] = T_k
    
    return K,F
    
            
            
v = np.ones(number_of_rectangle)*0.2 ; k_values = fill_in_k(v,65,0.2,p)
K = find_K(v)
F = find_F()
K,F = apply_boundary(K,F) ; T = np.linalg.solve(K,F)  

#################################################### Optimisation ###################################################


def get_K(nx,ny,v):
    """
    Finds stiffness matrix K.
    Args:
        nx,ny (int): elements in x- and y-direction. 
                    Total amount of elements is nx*ny, 
                    total amount of grid points is (nx+1)*(ny+1).
        v (np.ndarray): vector of length (nx*ny) containing metal fractions of grid elements. 
    Returns:
        K (np.ndarray): Stiffness Matrix
    """
    
    
    
    return find_K(v,nx)


def linear_system_solve(nx,ny,K,F):
    """
    Solves the system KT=F, taking into account the boundary conditions. 
    Args:
        nx,ny (int): elements in x- and y-direction. 
                    Total amount of elements is nx*ny, 
                    total amount of grid points is (nx+1)*(ny+1).
        K (np.ndarray): Stiffness matrix
        F (np.ndarray): Right-hand-side vector

    Returns: T (np.ndarray): Temperature vector as a solution to KT=F
    """
    
    return np.linalg.solve(K,F)  
    












# def f(v):
#     D = 0 
#     k_values = fill_in_k(v,k_max,k_min,p)
#     for e,rectangle in enumerate(rectangles):
#         k_e = k_values[e]  
#         T_loc = np.zeros((4,1))
#         for l in range(4):
#             T_loc[l][0]=T[int(rectangle[l])-1]
        
#         K = np.array(local_matrix)*k_e
#         D+= ((T_loc.T)@K@T_loc)[0][0]
          
#     return D



# def grad_f(v):
#     gradient_D = np.zeros(number_of_rectangle)
#     k_values = fill_in_k(v,k_max,k_min,p)
#     for e,rectangle in enumerate(rectangles):
#         k_e = k_values[e]  
#         T_loc = np.zeros((4,1))
#         for l in range(4):
#             T_loc[l][0]=T[int(rectangle[l])-1]
        
#         K = np.array(local_matrix)*k_e
#         gradient_D[e] = ((-1/2)*p*((v[e])**(p-1))*(k_max-k_min)*((T_loc.T)@K@T_loc)*(1/k_e))[0][0]
        
#     return gradient_D
    

# def constraint(v):
#     g = -(1/len(v))*np.sum(v)+0.4
#     return g

# bounds =[]
# for i in range(number_of_rectangle):
#     bounds.append([0,1])


# constraints = [{'type': 'ineq', 'fun': constraint}]
# result = minimize(f, v, jac=grad_f, bounds=bounds, constraints=constraints)


# print("Solution optimale :", result.x)
# print("Valeur de la fonction objectif à la solution :", result.fun)
# print("Succès de l'optimisation :", result.success)


#T_matrix = T.reshape((number_of_points, number_of_points))

# X_matrix = (np.array(result.x)).reshape((number_of_points-1, number_of_points-1))     
# plt.imshow(X_matrix, cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()

# plt.figure(figsize=(6, 5))
# plt.imshow(T_matrix, cmap='hot', origin='lower', extent=[0, L, 0, L])
# plt.colorbar(label="Température (K)")
# plt.title("Distribution de la température")
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# plt.show()


