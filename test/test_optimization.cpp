#include <iostream>
#include <Eigen/Sparse>
#include <cassert>
#include <iostream>
#include <random>

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
    std::random_device rd;
    std::mt19937 gen(rd());  // Generator with a random seed
    std::uniform_real_distribution<> dis(0.0, 1.0);
    Eigen::Matrix4d K0;
    K0 << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;
    int nx = 19;
    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(nx+1);
    Eigen::VectorXd x(nx*nx);
    
    std::vector<int> ps = {1,2,3,4,5};
    for (int ft=0;ft<3;ft++){
        for (int p=1;p<6;p++){
            for (int i = 0; i < nx*nx; ++i) {
                x(i) = dis(gen) * 0.5 + 0.5;  // Fill the initial state vector with random values between 0 and 1
            }
            x = (x.array() - x.mean()) + 0.4;  // Adjust the mean to 0.4
            x = x.cwiseMin(1.0).cwiseMax(0.0); // Ensure all values remain in [0,1]
            
            optimize(K0,x,0.4,nx,nx,p,rectangles,0.01,293,ft);
            // every element of x must be close to zero or one. 
            VectorXd x_rounded = x.array().round();
            std::cout << "ft="<<ft<<", p="<<p<<std::endl;
            std::cout << "maximum deviaton from 0 or 1: "<<(x_rounded - x).cwiseAbs().maxCoeff() << std::endl;
            std::cout << "mean deviation from 0 or 1:"<< (x_rounded - x).cwiseAbs().mean() << std::endl;
            //the mean deviation cannot be too high. 
            tf::compareTolerance((x_rounded - x).cwiseAbs().mean(),0.,0.2/p); 
        }
    }  
}

int main() {
    test_optimization_function();
    return 0;
}