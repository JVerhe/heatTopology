#include <Eigen/Sparse>
#include <iostream>
#include <csignal>
#include <cstdlib>
#include "optimization.hpp"
#include "meshHelper.hpp"
#include <random>
#include <chrono>

int main() {

    int p = 2;
    int output = 0;
    double  L = 0.01;
    double T_k = 293;
    double vol_frac = 0.4;
    Eigen::Matrix4d local_matrix;
    int max_repeats = 6;
    local_matrix << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;

    std::cout << "ft dim time" << std::endl;
    for (int ft = 0; ft <= 2; ft++) {
        for (int number_of_points = 10; number_of_points <= 100; number_of_points += 5) {
            double total_time = 0;
            for (int repeat = 0; repeat < max_repeats; repeat++) {
                int size = (number_of_points - 1) * (number_of_points - 1); // Grid size
                std::random_device rd;
                std::mt19937 gen(rd());  // Generator with a random seed
                std::uniform_real_distribution<> dis(0.0, 1.0);

                Eigen::VectorXd x(size);
                for (int i = 0; i < size; ++i) {
                    x(i) = dis(gen) * 0.5 + 0.5;  // Fill the initial state vector with random values between 0 and 1
                }
                x = (x.array() - x.mean()) + 0.4;  // Adjust the mean to 0.4
                x = x.cwiseMin(1.0).cwiseMax(0.0); // Ensure all values remain in [0,1]
                std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(number_of_points);

                std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                optimize(local_matrix, x, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, L, T_k, ft, output);
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                std::chrono::duration<double> diff = t2 - t1;
                total_time += diff.count();

            }
            std::cout << ft << " " << number_of_points << " " << total_time / max_repeats << std::endl;
        }
    }

    return 0;
}