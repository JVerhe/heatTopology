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
    int ft = 2;
    Eigen::Matrix4d K0;
    K0 << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;
    int nx = 19;
    int p = 4;
    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(nx+1);
    VectorXd x = VectorXd::Constant(nx * nx, 0.4);
    
    optimize(K0,x,0.4,nx,nx,p,rectangles,0.01,293,ft);

    // every element of x must be close to zero or one. 
    VectorXd x_rounded = x.array().round();
    //std::cout << "x=" << x << std::endl << "x_rounded = " << x_rounded << std::endl;
    std::cout << "maximum deviaton from 0 or 1: "<<(x_rounded - x).cwiseAbs().maxCoeff() << std::endl;
    std::cout << "mean deviation from 0 or 1 with penal="<<p<<":"<< (x_rounded - x).cwiseAbs().mean() << std::endl;
    //the mean deviation cannot be too high. 
    tf::compareTolerance((x_rounded - x).cwiseAbs().mean(),0.,0.2/p);   
}

int main() {
    test_optimization_function();
    return 0;
}