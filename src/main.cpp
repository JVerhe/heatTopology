#include <Eigen/Sparse>
#include <iostream>
#include "optimization.hpp"
#include "meshHelper.hpp"


int main() {
    double  L = 0.01;
    int p = 2;
    double T_k = 293;
    int number_of_points = 20;
    double vol_frac = 0.4;

    Eigen::Matrix4d local_matrix;
    local_matrix << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;


    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(number_of_points);
    int ft = 0;
    Eigen::VectorXd x = optimize(local_matrix, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, L, T_k, ft);

    save_result_to_file(x, "output/final.txt");

    return 0;
}