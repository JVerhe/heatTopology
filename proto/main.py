import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from math import ceil, sqrt
from helper_functions import *


L = 0.01 
outlet = 0.004    
number_of_points =10
number_of_rectangle = (number_of_points-1)*(number_of_points-1)
p = 3
k_min = 0.2
k_max = 65
T_k =  293

local_matrix = [[2/3,-1/6,-1/3,-1/6],[-1/6,2/3,-1/6,-1/3],[-1/3,-1/6,2/3,-1/6],[-1/6,-1/3,-1/6,2/3]]

rectangles = create_rectangle_and_mesh(number_of_points)
coordinates = create_coordinates(L,number_of_points)
boundary_points = filter_boundary_points_with_index(coordinates,L)


v = np.ones(number_of_rectangle)*0.2
k_values = fill_in_k(v,65,0.2,p)
K = find_K(v,rectangles,number_of_points,local_matrix,k_min,k_max,p)
F = find_F(rectangles,number_of_points,L)
K,F = apply_boundary(K,F,boundary_points,T_k) ; T = np.linalg.solve(K,F)  


#################################################### Optimisation ###################################################


x = optimization(
    K0=local_matrix,
    F=F,
    max_vol_frac=0.4,
    nx=number_of_points-1,
    ny=number_of_points-1,
    penal=p,
    rectangles=rectangles,
    L=L,
    boundary_temp=T_k,
    ft=0,
)















def grad_f(v):
    """
    This function accomplishes the same as adjointMethod and is therefore redundant.
    """
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


