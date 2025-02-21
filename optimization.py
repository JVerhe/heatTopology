import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
import matplotlib.pyplot as plt
from math import ceil, sqrt

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

    #Prepare filter
    iH = np.ones((nx*ny*(2*(ceil(rmin)-1)+1)**2,1))
    jH = np.ones(iH.shape())
    sH = np.zeros(iH.shape())
    H = coo_matrix(sH, (iH, jH))
    Hs =  np.sum(H,axis=1)
    #TODO: Prepare filter further

    #prepare for optimization iterations
    x = np.tile(max_vol_frac,(ny,nx))
    x.flatten(order='F') #needed? from matrix to vector
    xPhys = x
    loop = 0
    change = 1
    #start iteration
    while change > 0.01:
        loop += 1
        #FE-analysis
        K = find_K(xPhys,nx)
        U = np.linalg.solve(K,F)
        #Objective function and sensitivity analysis
        c = f(xPhys)
        dc = adjoint(T=U,v=xPhys,corners=corners,p=penal) #sensitivity of objective
        dv = np.ones((ny, nx))# Sensitivity of volume constraint
        dv.flatten(order='F')#from matrix to vector, is this okay?

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
        move = 0.2 #m in paper
        while (l2-l1)/(l1+l2) > 1e-3:
            lmid = 0.5*(l1+l2)
            xnew = np.maximum(0,np.maximum(x-move,np.minimum(x+move,x*sqrt(-dc/dv/lmid))))
            if ft == 1 or ft == 0:
                xPhys = xnew
            elif ft == 2:
                xPhys = (H@xnew.ravel(order='F'))/Hs
            if np.sum(xPhys.ravel(order='F')) > max_vol_frac*nx*ny:
                l1 = lmid
            else: l2 = lmid
        change = np.max(np.abs(xnew.ravel(order='F')-x.ravel(order='F')))
        x = xnew
        #print results
        print(f" It.: {loop:5d} Obj.: {c:11.4f} Vol.: {np.mean(xPhys):7.3f} ch.: {change:7.3f}")
        # Plot density distribution
        
        plt.imshow(1 - xPhys, cmap='gray', vmin=0, vmax=1)
        plt.axis('equal')
        plt.axis('off')
        plt.colorbar()  # Optional, if you want a color scale
        plt.show()

