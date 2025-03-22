import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/home/ari_prezo/Bureau/Project/team_02/proto')
from meshHelper import create_rectangle_and_mesh


def load_result_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        U = np.array([float(line.strip()) for line in lines])
    return U

with open("/home/ari_prezo/Bureau/Project/team_02/build/config/config.txt", 'r') as file:
    a = file.readline().strip()
    number_of_points = int(file.readline().strip())


plot = True

L = 0.01
h = (L)/(number_of_points-1)


T_computed = load_result_from_file("/home/ari_prezo/Bureau/Project/team_02/build/output/temperature_mms.txt")
T_computed_matrix = T_computed.reshape((number_of_points, number_of_points))

if plot==True :
    plt.figure(figsize=(6, 5))
    plt.imshow(T_computed_matrix, cmap='magma', origin='lower', extent=[0, L, 0, L])
    plt.colorbar(label="K")
    plt.title("Temperature Distribution")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.show()


T_verification = load_result_from_file("/home/ari_prezo/Bureau/Project/team_02/build/output/correct_temperature_mms.txt")
T_verification_matrix = T_verification.reshape((number_of_points, number_of_points))
Err  = T_verification_matrix - T_computed_matrix

if plot==True :
    plt.figure(figsize=(6, 5))
    plt.imshow(Err, cmap='gray_r', origin='lower', extent=[0, L, 0, L])
    plt.colorbar(label="Error (K)")
    plt.title("Error with the true solution")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.show()

xi_q = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
eta_q = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
w_q = np.ones(4)


def shape_functions(xi, eta):
    N1 = (1 - xi) * (1 - eta) / 4
    N2 = (1 + xi) * (1 - eta) / 4
    N3 = (1 + xi) * (1 + eta) / 4
    N4 = (1 - xi) * (1 + eta) / 4
    
    return np.array([N1, N2, N3, N4])


def compute_error():

    error_L2 = 0.0
    J = (h**2) / 9
      
    rectangles = create_rectangle_and_mesh(number_of_points)
    number_of_elements = np.sqrt(len(rectangles))
    
    for i in range(len(rectangles)):
        rectangle = rectangles[i]
        i1 = int(rectangle[0])
        i2 = int(rectangle[1])
        i3 = int(rectangle[2])
        i4 = int(rectangle[3])
        Uh_nodes = np.array([T_computed[i1],T_computed[i2],T_computed[i3],T_computed[i4]])
        U_nodes = np.array([T_verification[i1],T_verification[i2],T_verification[i3],T_verification[i4]])
        
        for xi in xi_q:
            for eta in eta_q:
                N = shape_functions(xi, eta)  
                U_q = np.dot(N, U_nodes) 
                Uh_q = np.dot(N, Uh_nodes) 
                Err_2 = (Uh_q-U_q)**2
                error_L2 += np.sum(Err_2)*J


    error_L2 = np.sqrt(error_L2)
    print("Error in L2 norm : " + str(error_L2))
    print("Ratio Error/h^2 : " +str(error_L2/(h**2)))
    print("Log of the error in L2 norm : " + str(np.log10(error_L2)))
    print("Log of h : " + str(np.log10(h)))

if __name__ == "__main__":
    
    print("Size of the mesh : "+ str(number_of_points*number_of_points))
    compute_error()
    
