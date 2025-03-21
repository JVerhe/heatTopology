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



double function_q(double L, double x, double y, double period) {
    double factor = (M_PI * period) / (L / 2);
    double sin_x = std::sin(factor * x);
    double sin_y = std::sin(factor * y);
    return 2 * (factor * factor) * sin_x * sin_y;
}

double function_q_with_k_variable(double L, double x, double y, double period) {
    double factor = (M_PI * period) / (L / 2);
    double sin_x = std::sin(factor * x);
    double cos_x = std::cos(factor * x);
    double sin_y = std::sin(factor * y);
    return 2 * (factor * factor) * sin_x * sin_y * x - factor * sin_y * cos_x;
}




double function_T(double L, double x, double y, double period) {
    double factor = (M_PI * period) / (L / 2);
    double sin_x = std::sin(factor * x);
    double sin_y = std::sin(factor * y);
    return sin_x * sin_y + 293;
}



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

        // Ajouter F_local aux indices correspondants
        for (int l = 0; l < 4; ++l) {
            int index_i = rectangles[count][l];
            F(index_i) += F_local;  // Ajout à l'indice correspondant dans F
        }
    }
    std::cout<< F.size() << std::endl;
    return F;
}



std::vector<Eigen::Vector3d> filter_boundary_points_with_index_mms(const Eigen::MatrixXd& coordinates, double L) {
    /**
     * Identifie tous les points situés sur les bords du domaine.
     * Chaque point de frontière est stocké sous forme d'un vecteur 3D (index, x, y).
     *
     * @param coordinates: Matrice (N x 2) contenant les coordonnées (x, y).
     * @param L: La longueur du domaine.
     * @return Un vecteur de Eigen::Vector3d où chaque entrée contient (index, x, y).
     */

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
 * @brief It's a classical solver of the system 
 * 
 * @param x the vector containing the initial solution. On return, it contains the solution. 
 * @param K0 The constant part of the local conductivity matrix.
 * @param max_vol_frac The maximum amount of volume percentage of metal on the plate.
 * @param nx The amount of elements in x-direction.
 * @param ny The amount of elements in y-direction.
 * @param penal The penalization factor in the SIMP method.
 * @param rectangles The numbering of neighbouring grid points for each element.
 * @param L The side length of the plate.
 * @param boundary_temp The fixed temperature at the outlets.
 * @param ft The filtering option: 0=no filtering, 1=sensitivity filtering, 2=density filtering.
 * 
 * @return void
 */
VectorXd solve_simple(
    const Eigen::MatrixXd& K0,
    Eigen::VectorXd& x,
    double max_vol_frac,
    int nx, int ny,
    double penal,
    const std::vector<std::vector<int>>& rectangles,
    double L,
    double boundary_temp,
    int ft
) {
    double E_min = 0.2;
    double E_0 = 65;
    double rmin = 0.04 * nx;

    int N_total_elements = nx * ny;
    int N_total_points = (nx + 1) * (ny + 1);
    int N_points_1D = nx + 1;

    assert(x.size() == nx * ny);
    VectorXd x_phys = x;

    double curr_obj = 1e6;
    int loop = 0;
    double change = 1;
    Eigen::MatrixXd coordinates = create_coordinates(L, N_points_1D);
    std::vector<Eigen::Vector3d> boundary_points = filter_boundary_points_with_index_mms(coordinates, L);
    
    int period = 1; 

    std::vector<double> q_rectangle = create_q_rectangle_middle(rectangles,coordinates,L,N_points_1D,period);
    SparseMatrix<double> K = find_K(x_phys, rectangles, N_points_1D, K0, 0.2, 65.0, penal);
    VectorXd F = find_F_mms(rectangles,N_points_1D,L,q_rectangle);

    auto result = apply_boundary_conditions_optimized(K, F, boundary_points, boundary_temp);
    K = result.first;
    F = result.second;

    VectorXd T = VectorXd::Zero(F.size());
    solve_sparse_lin_sys_LU(K, F, T);

    /// SAVE THE RESULTS ///
    Eigen::VectorXd temperature_value = Eigen::VectorXd::Map(T.data(), T.size());
    char filename_temperature[100] = "output/temperature_mms.txt";
    save_result_to_file_mms(temperature_value, filename_temperature);
    /////

    Eigen::VectorXd T_True = create_T_point(coordinates,N_points_1D,L,period);
    char filename_correct_temperature[100] = "output/correct_temperature_mms.txt";
    save_result_to_file_mms(T_True, filename_correct_temperature);
    Eigen::VectorXd Err = T-T_True;




    return T;
}