import math
import numpy as np

##### Generate the vector q(x,y) = -2*k*(pi/L)^2*(sin(pi/L*x)*(sin(pi/L*y))

## L is actually 0.05 because we work on half the plane
def function_q(L,x,y,period):
    factor = (math.pi*period)/(L/2)
    sin_x = math.sin(factor*x)
    sin_y = math.sin(factor*y)
    return 2*((factor)**2)*sin_x*sin_y

def function_q_with_k_variable(L,x,y,period):
    factor = (math.pi*period)/(L/2)
    sin_x = math.sin(factor*x)
    cos_x= math.cos(factor*x)
    sin_y = math.sin(factor*y)
    return 2*((factor)**2)*sin_x*sin_y*x - factor*sin_y*cos_x



def function_T(L,x,y,period):
    factor = (math.pi*period)/(L/2)
    sin_x = math.sin(factor*x)
    sin_y = math.sin(factor*y)
    return sin_x*sin_y + 293
    

## We suppose here that k is constant, 
#  so k is not an array but rather a double
def create_q_point(coordinates,L,number_of_points,period):
    q_point =np.zeros(number_of_points*number_of_points)
    for i in range(number_of_points):
        for j in range(number_of_points):
            index = i*number_of_points+j
            x = coordinates[index][0]
            y = coordinates[index][1]
            q_point[index] = function_q(L,x,y,period)
    return q_point


## We want to make on each rectangle the 
## mean of the q values on each corner
def create_q_rectangle(rectangles,q_point):
    q_rectangle = np.zeros(len(rectangles))
    ## We retrieve then q_values bottom left, bottom right, 
    ## up right and up left and we compute the mean
    count = 0 
    for rectangle in rectangles : 
        q_value_bl = q_point[int(rectangle[0])]
        q_value_br = q_point[int(rectangle[1])]
        q_value_ur = q_point[int(rectangle[2])]
        q_value_ul = q_point[int(rectangle[3])]
        q_rectangle[count] = (1/4)*(q_value_bl+q_value_br+q_value_ur+q_value_ul)
        count += 1 
    
    return q_rectangle


## We want to make on each rectangle the 
## mean of the q values on each corner
def create_q_rectangle_middle(rectangles,coordinates,L,number_of_point,period):
    q_rectangle = np.zeros(len(rectangles))
    h = (L/2)/(number_of_point-1)
    ## We retrieve then q_values bottom left, bottom right, 
    ## up right and up left and we compute the mean
    count = 0 
    for rectangle in rectangles : 
        index = int(rectangle[0])
        [x0 , y0] = coordinates[index]
        x0 += h/2
        y0 -= h/2
        q_rectangle[count] = function_q(L,x0,y0,period)
        count += 1 
    
    return q_rectangle


def create_q_rectangle_middle_k_variable(rectangles,coordinates,L,number_of_point,period):
    q_rectangle = np.zeros(len(rectangles))
    h = (L/2)/(number_of_point-1)
    ## We retrieve then q_values bottom left, bottom right, 
    ## up right and up left and we compute the mean
    count = 0 
    for rectangle in rectangles : 
        index = int(rectangle[0])
        [x0 , y0] = coordinates[index]
        x0 += h/2
        y0 -= h/2
        q_rectangle[count] = function_q_with_k_variable(L,x0,y0,period)
        count += 1 
    
    return q_rectangle
    

def find_F_mms(rectangles,number_of_points,L,q_rectangle):
    h = (L/2)/(number_of_points-1)
    F = np.zeros(number_of_points * number_of_points)
    count = 0 
    for rectangle in rectangles: 
        F_local = ((h**2)/(4))*q_rectangle[count]
        for l in range(4):
            index_i = int(rectangle[l])  
            F[index_i] += F_local 
        count+=1
    return F




def find_K_with_k_variable(v,rectangles,number_of_points,local_matrix,coordinates,k_min,k_max,penal):
    ############
    ## v : density vector with values between 0 and 1
    ## rectangles : list of the number of nodes for each rectangle
    ## local_matrix : the local matrix on one element (always the same)
    ## coordinates : the coordinates in one fourth of the chip
    ## number of points : the number of mesh point in one fourth of the mesh
    ## L : the length of the chip
    #####
    K = np.zeros((number_of_points*number_of_points,number_of_points*number_of_points))
    k_values = np.zeros(len(rectangles))
    h = 0.005/(number_of_points-1)
    for i in range(number_of_points-1):
        for j in range(number_of_points-1):
            index = i*(number_of_points-1) + j
            k_values[index] = j*h + h/2
    
    for e,rectangle in enumerate(rectangles):
        k_e = k_values[e]  
        for l in range(4):
            for m in range(4):
                index_i = int(rectangle[l])
                index_j = int(rectangle[m])
                value   = local_matrix[l][m]
                K[index_i,index_j] += k_e*value
                
    return K


## In this exemple we take all the boundary points
def filter_boundary_points_with_index_mms(coordinates, L):
    ############
    ## coordinates : the coordinates in one fourth of the chip
    ## number of points : the number of mesh point in one fourth of the mesh
    ## L : the length of the chip
    #####
    boundary_points = []
    for idx, point in enumerate(coordinates):
        x, y = point
        if ((x == 0 or x == L/2) or (y==0 or y==L/2)):
            boundary_points.append([idx, x, y])  # Stocke (index, x, y)

    return boundary_points
   

## We need to create the solution in every point to 
## compare with the numerical solution
def create_T_point(coordinates,k,L,number_of_points,period):
    T_point =np.zeros(number_of_points*number_of_points)
    for i in range(number_of_points):
        for j in range(number_of_points):
            index = i*number_of_points+j
            x = coordinates[index][0]
            y = coordinates[index][1]
            T_point[index] = function_T(L,x,y,period)
    return T_point
