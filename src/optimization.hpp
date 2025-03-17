#ifndef optimization_hpp
#define optimization_hpp

#include <Eigen/Sparse>
#include "meshSolver.hpp"
#include "optHelper.hpp"
#include <fstream>
#include <iostream>
#include <meshHelper.hpp>
#include <cassert>
#include <chrono>

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


double get_cummulative_sum(std::vector<double> durations, int final_idx) {
    double sum = 0.;
    for (int i = 0; i < final_idx; i++) {
        sum += durations[i];
    }
    return sum;
}

double get_delta(std::vector<double> v, int final_idx) {
    return v[0] - v[final_idx];
}


/**
 * @brief Performs a heuristic algorithm for heat topology optimization.
 *
 * This function heuristically optimizes the distribution of metal on a FEM-discretized plate.
 * It ensures elements to be either fully metal or fully plastic via the SIMP method.
 * Results are saved in folder 'output'.
 * Inspiration for this code is taken from:
 * Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011).
 * Efficient topology optimization in MATLAB using 88 lines of code.
 * Structural and Multidisciplinary Optimization, 43, 1-16.
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
 * @param visualize If visualize=0 -> benchmark the code, if visualize = 1 -> save intermediate results and call python script to plot
 *
 * @return void
 */
void optimize(
    const Eigen::MatrixXd& K0,
    Eigen::VectorXd& x,
    double max_vol_frac,
    int nx, int ny,
    double penal,
    const std::vector<std::vector<int>>& rectangles,
    double L,
    double boundary_temp,
    int ft,
    int visualize = 0
) {
    std::vector<double> objective_values;
    std::vector<double> temperature_values;
    double E_min = 0.2;
    double E_0 = 65;
    double rmin = 0.04 * nx;

    int N_total_elements = nx * ny;
    int N_total_points = (nx + 1) * (ny + 1);
    int N_points_1D = nx + 1;

    SparseMatrix<double> H;
    VectorXd Hs; //sum of rows of H
    if (ft != 0) create_sparse_matrix(nx, ny, rmin, H, Hs);

    assert(x.size() == nx * ny);
    VectorXd x_phys = x;

    double curr_obj = 1e6;
    int loop = 0;
    double change = 1;
    Eigen::MatrixXd coordinates = create_coordinates(L, N_points_1D);
    std::vector<Eigen::Vector3d> boundary_points = filter_boundary_points_with_index(coordinates, L);
    std::chrono::steady_clock::time_point t1;
    std::vector<double> loop_durations;

    if (visualize == 0) {
        t1 = std::chrono::steady_clock::now();
    }
    while (change > 0.01 && loop < 200) {
        std::chrono::steady_clock::time_point loop_start;
        if (visualize == 0) {
            loop_start = std::chrono::steady_clock::now();
        }
        loop++;

        SparseMatrix<double> K = find_K(x_phys, rectangles, N_points_1D, K0, 0.2, 65.0, penal);
        VectorXd F = find_F(rectangles, N_points_1D, L);

        auto result = apply_boundary(K, F, boundary_points, boundary_temp);
        K = result.first;
        F = result.second;
        VectorXd U = VectorXd::Zero(F.size());
        solve_sparse_lin_sys(K, F, U);

        double c = objective(x_phys, rectangles, U, K0, 0.2, 65, penal);
        VectorXd dc = VectorXd::Zero(x_phys.size());
        adjoint(dc, U, x_phys, rectangles, K0, penal);
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

        find_new_densities(nx, ny, x, x_phys, xnew, dc, dv, H, Hs, ft, max_vol_frac);

        change = (xnew - x).cwiseAbs().maxCoeff();
        x = xnew;

        if (visualize == 0) {
            std::chrono::steady_clock::time_point loop_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> loop_time = loop_end - loop_start;
            loop_durations.push_back(loop_time.count());
        }


        objective_values.push_back(c);
        temperature_values.push_back(U.maxCoeff());

        if (visualize == 1) {
            std::cout << "It.: " << loop << " Obj.: " << c << " Vol.: " << x_phys.mean() << " ch.: " << change << std::endl;
            if (loop >= 2 && (objective_values[loop - 1] < objective_values[loop - 2] * 0.9) || loop == 1) {
                char filename[100];
                sprintf(filename, "output/result_obj%.1f_iteration%d.txt", objective_values[loop - 1], loop);
                save_result_to_file(x, filename);
            }
        }
    }

    // TODO implement with loops instead of this atrocity
    if (visualize == 0) { // Print timing results to the console
        std::chrono::steady_clock::time_point tend = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = tend - t1;
        double time = diff.count();
        std::cout << "____________Results for optimization of grid size " << nx << " x " << ny << "____________" << std::endl;
        std::cout << "Finished optimization process after " << time << "s time and " << loop << " iterations" << std::endl;
        std::cout << "Average time per iteration step: " << time / loop << "s" << std::endl;
        std::cout << "Average iterations per second: " << loop / time << "/s" << std::endl;
        std::cout << std::endl;

        int frac5 = loop / 20;
        int frac10 = loop / 10;
        int frac25 = loop / 4;
        int frac50 = loop / 2;
        double elapsed5 = get_cummulative_sum(loop_durations, frac5);
        double elapsed10 = get_cummulative_sum(loop_durations, frac10);
        double elapsed25 = get_cummulative_sum(loop_durations, frac25);
        double elapsed50 = get_cummulative_sum(loop_durations, frac50);
        double elapsed100 = get_cummulative_sum(loop_durations, loop);
        std::cout << "Average time per loop: " << std::endl;
        std::cout << "First 5%: " << elapsed5 / frac5 << std::endl;
        std::cout << "First 10%: " << elapsed10 / frac10 << std::endl;
        std::cout << "First 25%: " << elapsed25 / frac25 << std::endl;
        std::cout << "First 50%: " << elapsed50 / frac50 << std::endl;
        std::cout << "First 100%: " << elapsed100 / loop << std::endl;
        std::cout << std::endl;
        double obj5 = get_delta(objective_values, frac5);
        double obj10 = get_delta(objective_values, frac10);
        double obj25 = get_delta(objective_values, frac25);
        double obj50 = get_delta(objective_values, frac50);
        double obj100 = get_delta(objective_values, loop);

        std::cout << "Delta(objective) / time" << std::endl;
        std::cout << "First 5%: " << obj5 / elapsed5 << "/s" << std::endl;
        std::cout << "First 10%: " << obj10 / elapsed10 << "/s" << std::endl;
        std::cout << "First 25%: " << obj25 / elapsed25 << "/s" << std::endl;
        std::cout << "First 50%: " << obj50 / elapsed50 << "/s" << std::endl;
        std::cout << "First 100%: " << obj100 / elapsed100 << "/s" << std::endl;
    }

    if (visualize == 1) {
        Eigen::VectorXd objective_value = Eigen::VectorXd::Map(objective_values.data(), objective_values.size());
        char filename_values[100] = "output/objective_values.txt";
        save_result_to_file(objective_value, filename_values);

        Eigen::VectorXd temperature_value = Eigen::VectorXd::Map(temperature_values.data(), temperature_values.size());
        char filename_temperature[100] = "output/temperature.txt";
        save_result_to_file(temperature_value, filename_temperature);

        char outputfile[100] = "output/density.txt";
        save_result_to_file(x, outputfile);
    }
}

#endif