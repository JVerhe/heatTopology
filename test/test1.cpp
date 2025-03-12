#include <iostream>
#include <Eigen/Sparse>

// Include the next 5 lines in your own files
#include <testFramework.hpp>
#include <../src/meshHelper.hpp>
#define check(cond) check(cond, __LINE__, __FILE__)
#define compare(val1, val2) compare(val1, val2, __LINE__, __FILE__)
#define compareTolerance(val1, val2, tol) compareTolerance(val1, val2, tol, __LINE__, __FILE__)

void exampleTest1() { // Should work
    tf::check(true);
    tf::check(1); // 1 is equal to true in c++
    tf::compare(1, 1); // Arguments should be of the exact same type
    tf::compareTolerance(1., 1.1, 0.2);

    Eigen::Vector3d v1;
    v1 << 1, 2, 3;
    Eigen::Vector3d v2;
    v2 << 1, 2, 3;

    tf::compareTolerance(0., (v1 - v2).norm(), 0.1);

    create_rectangle_and_mesh(50); // Functions of the project can still be called
}

int main() {
    exampleTest1();
    std::cout << "test1 ran succesfully" << std::endl;
    return 0;
}

