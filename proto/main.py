import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from meshHelper import *
from optimization import optimize


#################################################### Discretisation ###################################################

L = 0.01 
number_of_points = 15
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

f = open("results.txt", "w")
solution = ""
for value in x:
    solution += str(value) + " "
f.write(solution)
f.close()

