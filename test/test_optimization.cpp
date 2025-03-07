#include <iostream>
#include <Eigen/Sparse>
#include <cassert>

// Include the next lines in your own files
#include <testFramework.hpp>
#include <../src/meshHelper.hpp>
#include <../src/optimization.hpp>
#define check(cond) check(cond, __LINE__, __FILE__)
#define compare(val1, val2) compare(val1, val2, __LINE__, __FILE__)
#define compareTolerance(val1, val2, tol) compareTolerance(val1, val2, tol, __LINE__, __FILE__)

using namespace Eigen;

void test_optimization_function() 
{
    int ft = 0;
    Eigen::Matrix4d K0;
    K0 << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;
    int nx = 20;
    int p = 5;
    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(nx+1);
    VectorXd x = VectorXd::Constant(nx * nx, 0.4);
    
    optimize(x, K0,0.4,nx,nx,p,rectangles,0.01,293,ft);;

    // every element of x must be close to 0 or one. 
    VectorXd x_rounded = x.array().round();
    //std::cout << "x=" << x << std::endl << "x_rounded = " << x_rounded << std::endl;
    std::cout << "maximum deviaton from 0 or 1: "<<(x_rounded - x).cwiseAbs().maxCoeff() << std::endl;
    //the maximum deviation cannot be too high. 
    //tf::compareTolerance((x_rounded - x).cwiseAbs().maxCoeff(),0.,0.05);  
    //current configuration gets stuck in a loop: we need to counter this.
    // when are we happy with the SIMP penalization? 
}

int main() {
    test_optimization_function();
    return 0;
}