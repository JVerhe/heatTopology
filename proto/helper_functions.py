import numpy as np
from adjointMethod import adjoint

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


def optimization(K0,F,max_vol_frac,nx,ny,penal,rectangles,L,boundary_temp,ft):
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
         rectangles(nested list): 
            Nested list with first index the elements and second index the nodes of said element
        L (float): square side length
        boundary_temp (float): fixed temperature at certain boundaries
        ft (int): if 0: no filtering
                  if 1: sensitivity filtering
                  if 2: density filtering
       
    
    Returns: an optimal vector x containing the volume fractions
    """
    #Total elements is nx*ny, total grid points is (nx+1)*(ny+1)
    assert(nx==ny)
    #define problem-specific parameters
    E_min = 0.2
    E_0 = 65
    rmin = 0.04*nx
    
    # number_of_points = int(np.sqrt(nx) +1) why??
    N_total_elements = nx*ny
    N_total_points = (nx+1)*(ny+1)
    N_points_1D = nx+1

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
    x = np.ones(nx*ny)*max_vol_frac
    xPhys = x
    loop = 0
    change = 0.1
    coordinates = create_coordinates(L,N_points_1D)
    boundary_points=filter_boundary_points_with_index(coordinates,L)
    #start iteration
    while change > 0.01:
        loop += 1
        #FE-analysis
        K = find_K(xPhys,rectangles,N_points_1D,K0,E_min,E_0,penal)
        F = find_F(rectangles,N_points_1D,L)
        K,F = apply_boundary(K,F,boundary_points,boundary_temp)
        # print(K.shape) ; print(F.shape)
        U = np.linalg.solve(K,F)
        
        #Objective function and sensitivity analysis
        c = objective(xPhys,rectangles,U,K0,E_min,E_0,penal)
        dc = adjoint(U,xPhys,rectangles,K0,p=penal).ravel() #sensitivity of objective
        assert(np.shape(dc) == np.shape(xPhys))
        dv = np.ones(len(dc))# Sensitivity of volume constraint

        #Filtering/Modification of sensitivities
        """
        if ft==1: #sensitivity filtering
            dc = (H @ (x.ravel(order='F') * dc.ravel(order='F'))) / Hs / np.maximum(1e-3, x.ravel(order='F'))
        elif ft==2: #density filtering
            dc = H@(dc.ravel(order='F')/Hs)
            dv = H@(dv.ravel(order='F')/Hs)
        #else (ft==0): continue
        """

        #optimality criteria update of design variables and physical densities
        l1=0
        l2=1e9
        move = 0.1 #m in paper
        while (l2-l1)/(l1+l2) > 1e-3:
            lmid = 0.5*(l1+l2)
            #print("l2=",l2)
            B = np.sqrt(np.maximum(-dc/(dv*lmid),0)) #neg. values close to e_mach do not have np.sqrt()
            Bx = B*x
            xnew = np.zeros(len(x))
            for e in range(len(xnew)):
                if Bx[e] <= max(0,x[e]-move): xnew[e] = max(0,x[e]-move)
                elif Bx[e] >= min(1,x[e]-move): xnew[e] = min(1,x[e]+move)
                else: xnew[e] = Bx[e]

            xnew = np.maximum(0, np.where(np.maximum(x - move, 0) > x * B, np.maximum(x - move, 0), np.where(np.minimum(x + move, 1) < x * B, np.minimum(x + move, 1), x * B)))

            # xnew = np.maximum(0,np.maximum(x-move,np.minimum(x+move,x*B)))

            #print("npmin",np.minimum(x+move,x*B))
            xPhys = xnew
            if np.sum(xPhys) > max_vol_frac*nx*nx:
                #print("np.sum(xPhys) too big:",np.sum(xPhys))
                l1 = lmid 
            else: l2 = lmid
            """
            if ft == 1 or ft == 0:
                xPhys = xnew
            elif ft == 2:
                xPhys = (H@xnew.ravel(order='F'))/Hs
            if np.sum(xPhys) > max_vol_frac*(nx):
                l1 = lmid
            else: l2 = lmid
            """
        
        #print(xPhys)
        change = np.max(np.abs(xnew-x))
        x = xnew
        #print results
        print(f" It.: {loop} Obj.: {c} Vol.: {np.mean(xPhys)} ch.: {change}")
        # Plot density distribution
    return x