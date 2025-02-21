#include <iostream>
#include <Eigen/Dense>  // Include the dense matrix module

int main() {
    Eigen::Matrix2d mat;  // 2x2 matrix of doubles
    mat << 1, 2,
        3, 4;  // Initializing the matrix

    Eigen::Vector2d vec(5, 6);  // 2D vector

    Eigen::Vector2d result = mat * vec;  // Matrix-vector multiplication

    std::cout << "Result:\n" << result << std::endl;
    return 0;
}
