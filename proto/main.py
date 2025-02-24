import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.optimize import minimize
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix 
import matplotlib.pyplot as plt
from math import ceil, sqrt



L = 0.01 
outlet = 0.004    
number_of_points =10; h = L/(number_of_points-1)
number_of_rectangle = (number_of_points-1)*(number_of_points-1)
p = 2
k_min = 0.2
k_max = 65
T_k =  293

local_matrix = [[2/3,-1/6,-1/3,-1/6],[-1/6,2/3,-1/6,-1/3],[-1/3,-1/6,2/3,-1/6],[-1/6,-1/3,-1/6,2/3]]


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
                K[index_i,index_j] += k_e*value
                
    return K


def find_F():
    
    F = np.zeros(number_of_points * number_of_points)

    for rectangle in rectangles: 
        F_local = (h**2)*(10**7) / (4)
    
        for l in range(4):
            index_i = int(rectangle[l])  
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
K = find_K(v,number_of_points)
F = find_F()
K,F = apply_boundary(K,F) ; T = np.linalg.solve(K,F)  

#################################################### Optimisation ###################################################


def f(v):
    D = 0 
    k_values = fill_in_k(v,k_max,k_min,p)
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        T_loc = np.zeros((4,1))
        for l in range(4):
            T_loc[l][0]=T[int(rectangle[l])-1]
        
        K = np.array(local_matrix)*k_e
        D+= ((T_loc.T)@K@T_loc)[0][0]
          
    return D



def grad_f(v):
    gradient_D = np.zeros(number_of_rectangle)
    k_values = fill_in_k(v,k_max,k_min,p)
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        T_loc = np.zeros((4,1))
        for l in range(4):
            T_loc[l][0]=T[int(rectangle[l])]
        
        K = np.array(local_matrix)*k_e
        gradient_D[e] = ((-1/2)*p*((v[e])**(p-1))*(k_max-k_min)*((T_loc.T)@K@T_loc)*(1/k_e))[0][0]
        
    return gradient_D
    




def optimization(K0,F,max_vol_frac,nx,ny,penal,corners,ft):
    """
    Perform heat topology optimization with volume constraints via heuristic algorithm.
    Sensitivity or density filtering can be applied. 
    
    Args:
        K0 (np.ndarray): stiffness matrix of size 4x4
        F (np.ndarray): Right-hand-side vector for solving linear system. 
        max_vol_frac (float): maximum volume fraction metal/plastic
        nx (int): amount of elements in x-direction
        ny (int): amount of elements in y-direction
        penal (float): (exponent) penalization for constructing stiffness matrix.
                       Typically, penal=3
                        E_e = E_min + x_e^penal * (E_0 - E_min)
         corners(nested list): 
            Nested list with first index the elements and second index the nodes of said element
        ft (int): if 0: no filtering
                  if 1: sensitivity filtering
                  if 2: density filtering
       
    
    Returns: an optimal vector x containing the volume fractions
    """
    #Total elements is nx*ny, total grid points is (nx+1)*(ny+1)

    #define problem-specific parameters
    Emin = 0.2
    E0 = 65
    rmin = 0.04*nx
    
    number_of_points = int(np.sqrt(nx) +1)

    #Prepare filter
    #iH = np.ones((nx*ny*(2*(ceil(rmin)-1)+1)**2,1))
    #jH = np.ones(np.shape(iH))
    #sH = np.zeros(np.shape(iH))
    #H = sH 
    #H = coo_matrix(sH, (iH, jH))
    #Hs =  np.sum(H,axis=1)
    #TODO: Prepare filter further

    #prepare for optimization iterations
    #x = np.tile(max_vol_frac,(ny,nx))
    #x = x.flatten(order='F') #needed? from matrix to vector
    x = np.ones(nx)*0.2
    xPhys = x
    loop = 0
    change = 0.1
    #start iteration
    while change > 0.01:
        loop += 1
        #FE-analysis
        K = find_K(xPhys,number_of_points) ; F=find_F()
        K,F = apply_boundary(K,F)
        print(K.shape) ; print(F.shape)
        U = np.linalg.solve(K,F)
        
        #Objective function and sensitivity analysis
        c = f(xPhys)
        dc = grad_f(xPhys) #sensitivity of objective
        dv = np.ones(len(dc))# Sensitivity of volume constraint
        #dv = dv.flatten(order='F')#from matrix to vector, is this okay?

        #Filtering/Modification of sensitivities
        if ft==1:
            dc = (H @ (x.ravel(order='F') * dc.ravel(order='F'))) / Hs / np.maximum(1e-3, x.ravel(order='F'))
        elif ft==2:
            dc = H@(dc.ravel(order='F')/Hs)
            dv = H@(dv.ravel(order='F')/Hs)
        #else (ft==0): continue

        #optimality criteria update of design variables and physical densities
        l1=0
        l2=1e9
        move = 0.1 #m in paper
        count = 0 
        while (l2-l1)/(l1+l2) > 1e-3 and count < 1000:
            count = count +1
            lmid = 0.5*(l1+l2)
            print(l2)
            B = np.zeros(len(dc))
            for i in range(len(B)):
                if dv[i] !=0:
                    B[i] = np.sqrt(-dc[i]/(dv[i]*lmid))
                else : 
                    B[i] = 0
                
            xnew = np.maximum(0,np.maximum(x-move,np.minimum(x+move,x*B)))
            if ft == 1 or ft == 0:
                xPhys = xnew
            elif ft == 2:
                xPhys = (H@xnew.ravel(order='F'))/Hs
            if np.sum(xPhys) > max_vol_frac*(nx):
                l1 = lmid
            else: l2 = lmid
            
            print(np.sum(xPhys),max_vol_frac*(nx))
            print(xPhys)
        change = np.max(np.abs(xnew-x))
        x = xnew
        #print results
        print(f" It.: {loop:5d} Obj.: {c:11.4f} Vol.: {np.mean(xPhys):7.3f} ch.: {change:7.3f}")
        # Plot density distribution
    return x




x = optimization(local_matrix,F,0.4,number_of_rectangle,number_of_rectangle,p,rectangles,0)



# T_matrix = T.reshape((number_of_points, number_of_points))
# plt.figure(figsize=(6, 5))
# plt.imshow(T_matrix, cmap='hot', origin='lower', extent=[0, L, 0, L])
# plt.colorbar(label="Température (K)")
# plt.title("Distribution de la température")
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# plt.show()



# X_matrix = (np.array(result.x)).reshape((number_of_points-1, number_of_points-1))     
# plt.imshow(X_matrix, cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()



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


