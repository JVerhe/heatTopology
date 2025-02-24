import numpy as np
from scipy.sparse import coo_matrix, csr_matrix 
import matplotlib.pyplot as plt
from math import ceil, sqrt
from meshHelper import *
from optHelper import adjoint

def optimize(K0,F,max_vol_frac,nx,ny,penal,rectangles,L,boundary_temp,ft):
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
    change = 0.4
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