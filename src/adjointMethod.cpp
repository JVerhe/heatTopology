#include <iostream>
#include <Eigen/Dense>
#include <map>
#include "adjointMethod.hpp"

int main() {

    using namespace Eigen;
    std::map<int, std::vector<int>> corners;
    corners[0] = { 3, 4, 0, 1 };
    corners[1] = { 4, 5, 1, 2 };
    corners[2] = { 6, 7, 3, 4 };
    corners[3] = { 7, 8, 4, 5 };

    Vector<double, 9> T{ 273, 274, 271, 273, 272, 275, 271, 270, 272 };
    Vector<double, 4> v{ 0.5, 0.5, 0.5, 0.5 };
    adjoint(T, v, corners);

    return 0;
}
