import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from meshHelper import *
from optimization import optimize
import sys

# try:
#     number_of_points = int(sys.argv[1])
# except:
#     print("You must give a number of discretization points when executing this file: ex. python3 main.py 30")
#     exit()

#################################################### Discretisation ###################################################

number_of_points = 20
L = 0.01 
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

