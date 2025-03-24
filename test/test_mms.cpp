#include <iostream>
#include <../src/mms.cpp>
#include <Eigen/Sparse>

int main() {

    int number_of_points = 100;
    int p = 1;
    int ft = 1;

    double L = 0.01;
    double T_k = 293;
    double vol_frac = 0.4;

    int size = (number_of_points - 1) * (number_of_points - 1);
    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(number_of_points);

    Eigen::Matrix4d local_matrix;
    local_matrix << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;

    Eigen::VectorXd x(size);

    // Test for different values of k_constant
    mms(local_matrix, x, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, 2 * L, T_k, 0);
    mms(local_matrix, x, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, 2 * L, T_k, 1);
    return 0;
}