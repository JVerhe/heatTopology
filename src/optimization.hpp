#ifndef optimization_hpp
#define optimization_hpp

#include <Eigen/Sparse>
#include "meshSolver.hpp"
#include "optHelper.hpp"
#include <fstream>
#include <iostream>
#include <meshHelper.hpp>
#include <cassert>

using namespace Eigen;

/**
 * @brief Saves a given Eigen vector to a file.
 *
 * This function writes the contents of the vector `U` to a text file, 
 * storing one element per line. If the file cannot be opened, 
 * an error message is printed to standard error.
 *
 * @param U The Eigen vector to be saved.
 * @param filename The relative path of the file where the vector will be stored.
 * @return void
 */
void save_result_to_file(const Eigen::VectorXd& U, const std::string& filename) {
    std::ofstream file(filename);
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
 * @brief Performs a heuristic algorithm for heat topology optimization.
 * 
 * This function heuristically optimizes the distribution of metal on a FEM-discretized plate. 
 * It ensures elements to be either fully metal or fully plastic via the SIMP method. 
 * Inspiration for this code is taken from:  
 * Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011). 
 * Efficient topology optimization in MATLAB using 88 lines of code. 
 * Structural and Multidisciplinary Optimization, 43, 1-16.
 * 
 * @param K0 The constant part of the local conductivity matrix.
 * @param max_vol_frac The maximum amount of volume percentage of metal on the plate.
 * @param nx The amount of elements in x-direction. 
 * @param ny The amount of elements in y-direction. 
 * @param penal The penalization factor in the SIMP method. 
 * @param rectangles The numbering of neighbouring grid points for each element. 
 * @param L The side length of the plate. 
 * @param boundary_temp The fixed temperature at the outlets.
 * @param ft The filtering option: 0=no filtering, 1=sensitivity filtering, 2=density filtering. 
 */
Eigen::VectorXd optimize(
    const Eigen::MatrixXd& K0, 
    double max_vol_frac,
    int nx, int ny, 
    double penal, 
    const std::vector<std::vector<int>>& rectangles, 
    double L,
    double boundary_temp, 
    int ft 
){
    double E_min = 0.2;
    double E_0 = 65;
    double rmin = 0.04 * nx;

    int N_total_elements = nx * ny;
    int N_total_points = (nx + 1) * (ny + 1);
    int N_points_1D = nx + 1;

    SparseMatrix<double> H;
    VectorXd Hs; //sum of rows of H
    if (ft!=0) create_sparse_matrix(nx, ny, rmin, H, Hs);


    VectorXd x = VectorXd::Constant(nx * ny, max_vol_frac);
    VectorXd x_phys = x;
    VectorXd x_sol = x;

    double curr_obj = 1e6;
    int loop = 0;
    double change = 1;
    Eigen::MatrixXd coordinates = create_coordinates(L, N_points_1D);
    std::vector<Eigen::Vector3d> boundary_points = filter_boundary_points_with_index(coordinates, L);

    while (change > 0.01 && loop < 200) {
        loop++;

        Eigen::SparseMatrix<double> K = find_K(x_phys, rectangles, N_points_1D, K0, 0.2, 65.0, penal);
        Eigen::VectorXd F = find_F(rectangles, N_points_1D, L);

        auto result = apply_boundary(K, F, boundary_points, boundary_temp);
        K = result.first;
        F = result.second;
        VectorXd U = VectorXd::Zero(F.size());
        solve_sparse_lin_sys(K,F,U);

        double c = objective(x_phys, rectangles, U, K0, 0.2, 65, penal);
        Eigen::VectorXd dc = adjoint(U, x_phys, rectangles, K0, penal);
        VectorXd dv = VectorXd::Ones(dc.size());

        if (ft == 1) { // Sensitivity filtering
            dc = H * (x.array() * dc.array()).matrix();
            dc = dc.array() / (Hs.array().max(1e-3));
        }
        else if (ft == 2) { // Density filtering
            dc = H * (dc.array() / Hs.array()).matrix();
            dv = H * (dv.array() / Hs.array()).matrix();
        }
        Eigen::VectorXd xnew = Eigen::VectorXd::Zero(x.size());

        find_new_densities(nx,ny,x,x_phys,xnew,dc,dv,H,Hs,ft,max_vol_frac);

        change = (xnew - x).cwiseAbs().maxCoeff();
        x = xnew;

        if (c < curr_obj){ //save solution if it is the best one yet
            curr_obj = c;
            x_sol = x;
        }

        std::cout << "It.: " << loop << " Obj.: " << c << " Vol.: " << x_phys.mean() << " ch.: " << change << std::endl;
        if (loop % 20 == 0) {
            char filename[100];
            sprintf(filename, "output/result_p%d_itteration%d.txt", penal, loop);
            save_result_to_file(x, filename);
        }
    }
    return x_sol;
}
#endif