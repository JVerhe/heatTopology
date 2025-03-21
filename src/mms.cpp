#include <Eigen/Sparse>
#include "meshSolver.hpp"
#include "optHelper.hpp"
#include <fstream>
#include <iostream>
#include <meshHelper.hpp>
#include <cassert>
#include <iomanip>  
#include <limits>  


using namespace Eigen;



void save_result_to_file_mms(const Eigen::VectorXd& U, const std::string& filename) {
    std::ofstream file(filename);
    file << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    if (file.is_open()) {
        for (int i = 0; i < U.size(); ++i) {
            file << U(i) << std::endl;
        }
        file.close();
        std::cout << "Saved to build/" << filename << std::endl;
    }
    else {
        std::cerr << "Error when saving " << filename << std::endl;
    }
}


/**
 * @brief Computes the heat source term q at a given point (x, y).
 * 
 * This function calculates the value of the heat source term q using a sinusoidal 
 * pattern based on the given length and period. It assumes that the thermal 
 * conductivity k is uniformly equal to 1 throughout the domain.
 * 
 * @param L The characteristic length of the domain.
 * @param x The x-coordinate where the function is evaluated.
 * @param y The y-coordinate where the function is evaluated.
 * @param period The number of oscillations in the domain.
 * 
 * @return The computed value of the heat source term at (x, y).
 */
double function_q(double L, double x, double y, double period) {
    double factor = (M_PI * period) / (L / 2);
    double sin_x = std::sin(factor * x);
    double sin_y = std::sin(factor * y);
    return 2 * (factor * factor) * sin_x * sin_y;
}

/**
 * @brief Computes the heat source term q at a given point (x, y) for a spatially varying conductivity.
 * 
 * This function calculates the value of the heat source term q using a sinusoidal 
 * pattern based on the given length and period. It assumes that the thermal 
 * conductivity k varies spatially and is given by k(x, y) = x.
 * 
 * @param L The characteristic length of the domain.
 * @param x The x-coordinate where the function is evaluated.
 * @param y The y-coordinate where the function is evaluated.
 * @param period The number of oscillations in the domain.
 * 
 * @return The computed value of the heat source term at (x, y).
 */
double function_q_with_k_variable(double L, double x, double y, double period) {
    double factor = (M_PI * period) / (L / 2);
    double sin_x = std::sin(factor * x);
    double cos_x = std::cos(factor * x);
    double sin_y = std::sin(factor * y);
    return 2 * (factor * factor) * sin_x * sin_y * x - factor * sin_y * cos_x;
}


/**
 * @brief Computes the manufactured solution T at a given point (x, y).
 * 
 * This function calculates the manufactured solution T using a sinusoidal 
 * pattern based on the given length and period. The solution is designed 
 * to satisfy a predefined heat equation and includes a baseline temperature of 293 K.  
 * 
 * In this case, we choose to impose homogeneous boundary conditions everywhere.
 * 
 * @param L The characteristic length of the domain.
 * @param x The x-coordinate where the function is evaluated.
 * @param y The y-coordinate where the function is evaluated.
 * @param period The number of oscillations in the domain.
 * 
 * @return The computed temperature T at (x, y).
 */
double function_T(double L, double x, double y, double period) {
    double factor = (M_PI * period) / (L / 2);
    double sin_x = std::sin(factor * x);
    double sin_y = std::sin(factor * y);
    return sin_x * sin_y + 293;
}

/**
 * @brief Computes the heat source term q at the center of each rectangle for k = 1.
 * 
 * This function evaluates the heat source term q at the midpoint of each rectangle 
 * in the computational grid. The midpoint is determined based on the coordinates 
 * of one of the rectangle's corners, adjusted by half the grid spacing h. 
 * It assumes that the thermal conductivity k is uniformly equal to 1.
 * 
 * @param rectangles A list of rectangles, where each entry contains the indices of its corner points.
 * @param coordinates A matrix containing the (x, y) coordinates of all points in the grid.
 * @param L The characteristic length of the domain.
 * @param number_of_point The number of points along one dimension of the grid.
 * @param period The number of oscillations in the domain.
 * 
 * @return A vector containing the computed values of q at the center of each rectangle.
 */
std::vector<double> create_q_rectangle_middle(
    const std::vector<std::vector<int>>& rectangles,
    const Eigen::MatrixXd& coordinates,
    double L,
    int number_of_point,
    double period) 
{
    std::vector<double> q_rectangle(rectangles.size(), 0.0);
    double h = (L / 2) / (number_of_point - 1);

    for (size_t count = 0; count < rectangles.size(); ++count) {
        int index = rectangles[count][0];
        double x0 = coordinates(index, 0) + h / 2;
        double y0 = coordinates(index, 1) - h / 2;
        q_rectangle[count] = function_q(L, x0, y0, period);
    }
    return q_rectangle;
}


/**
 * @brief Computes the heat source term q at the center of each rectangle for k(x, y) = x.
 * 
 * This function evaluates the heat source term q at the midpoint of each rectangle 
 * in the computational grid. The midpoint is determined based on the coordinates 
 * of one of the rectangle's corners, adjusted by half the grid spacing h.  
 * It assumes that the thermal conductivity k varies spatially and is given by k(x, y) = x.
 * 
 * @param rectangles A list of rectangles, where each entry contains the indices of its corner points.
 * @param coordinates A matrix containing the (x, y) coordinates of all points in the grid.
 * @param L The characteristic length of the domain.
 * @param number_of_point The number of points along one dimension of the grid.
 * @param period The number of oscillations in the domain.
 * 
 * @return A vector containing the computed values of q at the center of each rectangle.
 */
std::vector<double> create_q_rectangle_middle_k_variable(
    const std::vector<std::vector<int>>& rectangles,
    const Eigen::MatrixXd& coordinates,
    double L,
    int number_of_point,
    double period) 
{
    std::vector<double> q_rectangle(rectangles.size(), 0.0);
    double h = (L / 2) / (number_of_point - 1);

    for (size_t count = 0; count < rectangles.size(); ++count) {
        int index = rectangles[count][0];
        double x0 = coordinates(index, 0) + h / 2;
        double y0 = coordinates(index, 1) - h / 2;
        q_rectangle[count] = function_q_with_k_variable(L, x0, y0, period);
    }
    return q_rectangle;
}


/**
 * @brief Computes the stiffness matrix with modified thermal conductivity k(x, y) = x.
 * 
 * This function computes the stiffness matrix for a system, where the thermal 
 * conductivity k is spatially variable and is given by \( k(x, y) = x \). 
 * The conductivity values are calculated using the grid coordinates and the 
 * density vector \( v \). The function then builds the stiffness matrix based 
 * on the local element matrices and these computed conductivity values.
 * 
 * @param v The density vector used for calculating the conductivity values.
 * @param rectangles A list of rectangles, where each entry contains the indices 
 *        of its corner points in the global system.
 * @param number_of_points The number of points along one dimension of the grid.
 * @param local_matrix The local stiffness matrix for each rectangle.
 * @param k_min The minimum value for the conductivity (not directly used here, but 
 *        may be helpful for boundary conditions or scaling).
 * @param k_max The maximum value for the conductivity (not directly used here, but 
 *        may be helpful for boundary conditions or scaling).
 * @param penal The penalization factor (used in the SIMP method for controlling 
 *        material distribution, not directly used in this function).
 * @param L The characteristic length of the domain.
 * 
 * @return A sparse matrix representing the global stiffness matrix K.
 */
Eigen::SparseMatrix<double> find_K_with_k_variable(
    const Eigen::VectorXd& v,
    const std::vector<std::vector<int>>& rectangles,
    int number_of_points,
    const Eigen::Matrix4d& local_matrix,
    double k_min,
    double k_max,
    double penal,
    double L
) {
    using T = Eigen::Triplet<double>;
    std::vector<T> triplets;

    Eigen::VectorXd k_values(rectangles.size());
    double h = (L/2) / (number_of_points - 1);
    for (int i = 0; i < number_of_points - 1; ++i) {
        for (int j = 0; j < number_of_points - 1; ++j) {
            int index = i * (number_of_points - 1) + j;
            k_values[index] = j * h + h / 2;
        }
    }

    int size = number_of_points * number_of_points;
    Eigen::SparseMatrix<double> K(size, size);


    for (size_t e = 0; e < rectangles.size(); ++e) {
        double k_e = k_values[e];

        for (int l = 0; l < 4; ++l) {
            for (int m = 0; m < 4; ++m) {
                int index_i = rectangles[e][l];
                int index_j = rectangles[e][m];
                double value = local_matrix(l, m);

                triplets.emplace_back(index_i, index_j, k_e * value);
            }
        }
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}


/**
 * @brief Computes the nodal force vector for a given set of rectangles and heat source term.
 * 
 * This function calculates the global nodal force vector based on the heat source term 
 * applied to each rectangle element. It can handle both constant and variable heat source 
 * terms, depending on the values provided in the `q_rectangle` vector. The force vector is 
 * computed by integrating the heat source term over the area of each element, which is assumed 
 * to have a rectangular shape.
 * 
 * The force is calculated using the formula \( F_{\text{local}} = \frac{h^2}{4} \cdot q_{\text{rectangle}}[i] \), 
 * where \( h \) is the grid spacing and \( q_{\text{rectangle}}[i] \) is the heat source term at 
 * the center of the corresponding rectangle.
 * 
 * @param rectangles A list of rectangles, where each entry contains the indices of its corner points.
 * @param number_of_points The number of points along one dimension of the grid.
 * @param L The characteristic length of the domain.
 * @param q_rectangle A vector containing the heat source term evaluated at the center of each rectangle.
 * 
 * @return A vector containing the global nodal forces, with one entry per node.
 */
Eigen::VectorXd find_F_mms(
    const std::vector<std::vector<int>>& rectangles, 
    int number_of_points, 
    double L, 
    const std::vector<double>& q_rectangle) 
{
    double h = (L / 2) / (number_of_points - 1);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(number_of_points * number_of_points);  // Initialisation à zéro

    for (size_t count = 0; count < rectangles.size(); ++count) {
        double F_local = ((h * h) / 4.0) * q_rectangle[count];

        for (int l = 0; l < 4; ++l) {
            int index_i = rectangles[count][l];
            F(index_i) += F_local;  
        }
    }
    std::cout<< F.size() << std::endl;
    return F;
}


/**
 * @brief Selects the boundary points for Dirichlet conditions in the MMS case.
 * 
 * This function selects the boundary points for Dirichlet conditions, which are 
 * defined as all points located on the edges of the square domain. In the context 
 * of the Method of Manufactured Solutions (MMS), these boundary points are treated 
 * as fixed, and Dirichlet conditions are applied at these points.
 * 
 * Specifically, the boundary points are selected as those where the \( x \)- or 
 * \( y \)-coordinates are either 0 or \( L/2 \), where \( L/2 \) is the length of 
 * the domain. These boundary points are then returned in a vector of \( \text{Eigen::Vector3d} \), 
 * with each vector containing the index and coordinates of the boundary point.
 * 
 * @param coordinates The matrix of all points' coordinates in the grid, with each row representing 
 *        a point's \( (x, y) \) coordinates.
 * @param L The characteristic length of the square domain.
 * 
 * @return A vector containing the boundary points, each represented by its index and coordinates.
 */
std::vector<Eigen::Vector3d> filter_boundary_points_with_index_mms(const Eigen::MatrixXd& coordinates, double L) {

    std::vector<Eigen::Vector3d> boundary_points;

    for (int idx = 0; idx < coordinates.rows(); ++idx) {
        double x = coordinates(idx, 0);
        double y = coordinates(idx, 1);

        if (x == 0 || x == L / 2 || y == 0 || y == L / 2) {
            boundary_points.emplace_back(idx, x, y);
        }
    }

    return boundary_points;
}

/**
 * @brief Evaluates the exact solution \( T \) at each point of the mesh.
 * 
 * This function computes the exact temperature field \( T \) at each point of the mesh 
 * based on the manufactured solution. The function evaluates \( T(x, y) \) at each grid point 
 * and returns the results in a vector. The exact solution is assumed to be of the form \( T(x, y) = 
 * \sin\left(\frac{\pi \cdot x}{L}\right) \sin\left(\frac{\pi \cdot y}{L}\right) + 293 \), where \( L \) 
 * is the domain length and the period is used for the sinusoidal components.
 * 
 * @param coordinates The matrix (N x 2) containing the coordinates \( (x, y) \) for each point in the mesh.
 * @param number_of_points The number of points in each direction of the mesh.
 * @param L The characteristic length of the domain.
 * @param period The period of the sinusoidal signal used in the exact solution.
 * 
 * @return A vector containing the values of \( T \) at each point in the mesh.
 */
Eigen::VectorXd create_T_point(const Eigen::MatrixXd& coordinates, int number_of_points, double L, double period) {
    /**
     * Génère la solution T en chaque point du maillage pour la comparer avec la solution numérique.
     *
     * @param coordinates: Matrice (N x 2) contenant les coordonnées (x, y).
     * @param number_of_points: Nombre de points dans chaque direction.
     * @param L: Longueur du domaine.
     * @param period: Période du signal.
     * @return Un vecteur contenant les valeurs de T en chaque point du maillage.
     */

    Eigen::VectorXd T_point(number_of_points * number_of_points);

    for (int i = 0; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points; ++j) {
            int index = i * number_of_points + j;
            double x = coordinates(index, 0);
            double y = coordinates(index, 1);
            T_point[index] = function_T(L, x, y, period);
        }
    }

    return T_point;
}


/**
 * @brief Performs a verification of the code using the Method of Manufactured Solutions (MMS).
 * 
 * This function verifies the correctness of the code by evaluating the exact solution \( T \) 
 * at each point of the mesh, and comparing it with the computed solution. The function is 
 * designed to test the implementation of the finite element method with either a constant 
 * conductivity \( k = 1 \) or a spatially varying conductivity \( k(x, y) = x \). The 
 * `k_constant` parameter is used to specify which case to use: if `k_constant` is true, 
 * the conductivity is assumed to be constant with a value of 1, otherwise, it takes the 
 * form \( k(x, y) = x \).
 * 
 * This allows for a controlled comparison of the numerical solution with the exact solution 
 * in the context of MMS, ensuring that the code produces correct results for both types 
 * of conductivity.
 * 
 * @param coordinates The matrix (N x 2) containing the coordinates \( (x, y) \) for each point in the mesh.
 * @param number_of_points The number of points in each direction of the mesh.
 * @param L The characteristic length of the domain.
 * @param period The period of the sinusoidal signal used in the exact solution.
 * @param k_constant A boolean that determines the type of conductivity: 
 *        `true` for constant conductivity (k = 1), 
 *        `false` for spatially varying conductivity (k(x, y) = x).
 * 
 * @return A vector containing the values of \( T \) at each point in the mesh.
 */
void mms(
    const Eigen::MatrixXd& K0,
    Eigen::VectorXd& x,
    double max_vol_frac,
    int nx, int ny,
    double penal,
    const std::vector<std::vector<int>>& rectangles,
    double L,
    double boundary_temp,
    bool k_constant

) {
    double E_min = 0.2;
    double E_0 = 65;
    double rmin = 0.04 * nx;

    int N_total_elements = nx * ny;
    int N_total_points = (nx + 1) * (ny + 1);
    int N_points_1D = nx + 1;

    for (int i = 0; i < rectangles.size() ; ++i) {x(i) = 0.8/(65-0.2);  }
    VectorXd x_phys = x;

    double curr_obj = 1e6;
    int loop = 0;
    double change = 1;
    Eigen::MatrixXd coordinates = create_coordinates(L, N_points_1D);
    std::vector<Eigen::Vector3d> boundary_points = filter_boundary_points_with_index_mms(coordinates, L);
    
    int period = 1; 

    SparseMatrix<double> K;
    VectorXd F;

    if(k_constant){
        std::vector<double> q_rectangle = create_q_rectangle_middle(rectangles,coordinates,L,N_points_1D,period);
        K = find_K(x_phys, rectangles, N_points_1D, K0, 0.2, 65.0, penal);
        F = find_F_mms(rectangles,N_points_1D,L,q_rectangle);
    }else{
        std::vector<double> q_rectangle = create_q_rectangle_middle_k_variable(rectangles,coordinates,L,N_points_1D,period);
        K = find_K_with_k_variable(x_phys, rectangles, N_points_1D, K0, 0.2, 65.0, penal,L);
        F = find_F_mms(rectangles,N_points_1D,L,q_rectangle);
    }

    auto result = apply_boundary_conditions_optimized(K, F, boundary_points, boundary_temp);
    K = result.first;
    F = result.second;

    VectorXd T = VectorXd::Zero(F.size());
    solve_sparse_lin_sys_LU(K, F, T);

    ////////////////////// SAVE THE RESULTS //////////////////////////////////////////
    char filename_temperature[100] = "output/temperature_mms.txt";                  //
    save_result_to_file_mms(T, filename_temperature);                               //
                                                                                    // 
    Eigen::VectorXd T_True = create_T_point(coordinates,N_points_1D,L,period);      //
    char filename_correct_temperature[100] = "output/correct_temperature_mms.txt";  //
    save_result_to_file_mms(T_True, filename_correct_temperature);                  //
    //////////////////////////////////////////////////////////////////////////////////

    return ;
}