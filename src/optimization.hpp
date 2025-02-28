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


Eigen::VectorXd optimize(
    const Eigen::MatrixXd& K0, 
    double max_vol_frac,
    int nx, int ny, 
    double penal, 
    const std::vector<std::vector<int>>& rectangles, 
    double L,
    double boundary_temp, 
    int ft) {

    double E_min = 0.2;
    double E_0 = 65;
    double rmin = 0.04 * nx;

    int N_total_elements = nx * ny;
    int N_total_points = (nx + 1) * (ny + 1);
    int N_points_1D = nx + 1;

    SparseMatrix<double> H;
    VectorXd Hs; //sum of rows of H
    if (ft!=0) createSparseMatrix(nx, ny, rmin, H, Hs);


    VectorXd x = VectorXd::Constant(nx * ny, max_vol_frac);
    VectorXd xPhys = x;
    int loop = 0;
    double change = 1;
    Eigen::MatrixXd coordinates = create_coordinates(L, N_points_1D);
    std::vector<Eigen::Vector3d> boundary_points = filter_boundary_points_with_index(coordinates, L);

    while (change > 0.01 && loop < 200) {
        loop++;

        Eigen::SparseMatrix<double> K = find_K(xPhys, rectangles, N_points_1D, K0, 0.2, 65.0, penal);
        Eigen::VectorXd F = find_F(rectangles, N_points_1D, L);

        auto result = apply_boundary(K, F, boundary_points, boundary_temp);
        K = result.first;
        F = result.second;
        VectorXd U = VectorXd::Zero(F.size());
        solve_sparse_lin_sys(K,F,U);

        double c = objective(xPhys, rectangles, U, K0, 0.2, 65, penal);
        Eigen::VectorXd dc = adjoint(U, xPhys, rectangles, K0, penal);
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

        find_new_densities(nx,ny,x,xPhys,xnew,dc,dv,H,Hs,ft,max_vol_frac);

        change = (xnew - x).cwiseAbs().maxCoeff();
        x = xnew;

        std::cout << "It.: " << loop << " Obj.: " << c << " Vol.: " << xPhys.mean() << " ch.: " << change << std::endl;
        if (loop % 20 == 0) {
            char filename[100];
            sprintf(filename, "output/result_p%d_itteration%d.txt", penal, loop);
            save_result_to_file(x, filename);
        }
    }
    return x;
}

#endif