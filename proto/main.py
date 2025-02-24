import numpy as np
import matplotlib.pyplot as plt
from meshHelper import *
from optimization import optimize


L = 0.01 
number_of_points = 9
p = 3
T_k =  293

local_matrix = [[2/3,-1/6,-1/3,-1/6],[-1/6,2/3,-1/6,-1/3],[-1/3,-1/6,2/3,-1/6],[-1/6,-1/3,-1/6,2/3]]

rectangles = create_rectangle_and_mesh(number_of_points)
coordinates = create_coordinates(L,number_of_points)
boundary_points = filter_boundary_points_with_index(coordinates,L)
F = find_F(rectangles,number_of_points,L)


#################################################### Optimisation ###################################################

x = optimize(
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

print(x)
for e in x: assert(e<=1)


# T_matrix = T.reshape((number_of_points, number_of_points))
# plt.figure(figsize=(6, 5))
# plt.imshow(T_matrix, cmap='hot', origin='lower', extent=[0, L, 0, L])
# plt.colorbar(label="Température (K)")
# plt.title("Distribution de la température")
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# plt.show()



X_matrix = (np.array(x)).reshape((number_of_points-1, number_of_points-1))     
plt.imshow(X_matrix, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()



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


