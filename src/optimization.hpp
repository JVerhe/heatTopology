#ifndef optimization_hpp
#define optimization_hpp

#include <Eigen/Sparse>
#include "meshSolver.hpp"
#include "optHelper.hpp"
#include <fstream>
#include <iostream>
#include <meshHelper.hpp>

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

Eigen::VectorXd optimize(const Eigen::MatrixXd& K0, double max_vol_frac,
    int nx, int ny, double penal, const std::vector<std::vector<int>>& rectangles, double L,
    double boundary_temp, int ft) {

    double E_min = 0.2;
    double E_0 = 65;
    double rmin = 0.04 * nx;

    int N_total_elements = nx * ny;
    int N_total_points = (nx + 1) * (ny + 1);
    int N_points_1D = nx + 1;

    std::vector<int> iH, jH;
    std::vector<double> sH;

    for (int i1 = 0; i1 < nx; ++i1) {
        for (int j1 = 0; j1 < ny; ++j1) {
            int e1 = i1 * ny + j1;
            for (int i2 = std::max(i1 - static_cast<int>(std::ceil(rmin)) + 1, 0); i2 < std::min(i1 + static_cast<int>(std::ceil(rmin)), nx); ++i2) {
                for (int j2 = std::max(j1 - static_cast<int>(std::ceil(rmin)) + 1, 0); j2 < std::min(j1 + static_cast<int>(std::ceil(rmin)), ny); ++j2) {
                    int e2 = i2 * ny + j2;
                    double weight = std::max(0.0, rmin - std::sqrt((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2)));
                    iH.push_back(e1);
                    jH.push_back(e2);
                    sH.push_back(weight);
                }
            }
        }
    }

    SparseMatrix<double> H;
    VectorXd Hs;
    createSparseMatrix(nx, ny, iH, jH, sH, H, Hs);


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


        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
        solver.compute(K);

        solver.setMaxIterations(1000);
        solver.setTolerance(1e-4);

        Eigen::VectorXd U = solver.solve(F);

        double c = objective(xPhys, rectangles, U, K0, 0.2, 65, penal);
        Eigen::VectorXd dc = adjoint(U, xPhys, rectangles, K0, penal);

        VectorXd dv = VectorXd::Ones(dc.size());

        if (ft == 1) {
            dc = H * (x.array() * dc.array()).matrix();
            dc = dc.array() / (Hs.array().max(1e-3));
        }
        else if (ft == 2) { // Density filtering
            dc = H * (dc.array() / Hs.array()).matrix();
            dv = H * (dv.array() / Hs.array()).matrix();
        }


        double l1 = 0;
        double l2 = 1e9;
        double move = 0.4;
        Eigen::VectorXd xnew = Eigen::VectorXd::Zero(x.size());
        while ((l2 - l1) / (l2 + l1) > 0.001) {

            double lmid = 0.5 * (l1 + l2);
            VectorXd B = (-dc.array() / (dv.array() * lmid)).max(0).sqrt();
            VectorXd Bx = B.array() * x.array();




            xnew = Eigen::VectorXd::Zero(x.size());
            for (int e = 0; e < x.size(); ++e) {
                if (Bx[e] <= std::max(0.0, x[e] - move)) {
                    xnew[e] = std::max(0.0, x[e] - move);
                }
                else if (Bx[e] >= std::min(1.0, x[e] + move)) {
                    xnew[e] = std::min(1.0, x[e] + move);
                }
                else {
                    xnew[e] = Bx[e];
                }
            }


            for (size_t i = 0; i < x.size(); ++i) {
                double x_val = x[i];
                double B_val = B[i];

                // Calcul de la valeur de xnew[i]
                double max_x_move = std::max(x_val - move, 0.0);
                double min_x_move = std::min(x_val + move, 1.0);
                double x_times_B = x_val * B_val;

                if (max_x_move > x_times_B) {
                    xnew[i] = std::max(max_x_move, 0.0);
                }
                else if (min_x_move < x_times_B) {
                    xnew[i] = std::min(min_x_move, 1.0);
                }
                else {
                    xnew[i] = x_times_B;
                }

                xnew[i] = std::max(xnew[i], 0.0);
            }


            if (ft == 0 || ft == 1) {
                xPhys = xnew;
            }
            else if (ft == 2) {
                xPhys = H * xnew;
                xPhys = xPhys.cwiseQuotient(Hs);
            }

            if (xPhys.sum() > max_vol_frac * nx * nx) {
                l1 = lmid;
            }
            else {
                l2 = lmid;
            }
        }
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