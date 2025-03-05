#include <iostream>
#include <Eigen/Dense>
#include <cassert>

// Include the next 5 lines in your own files
#include <testFramework.hpp>
#include <../src/optHelper.hpp>
#define check(cond) check(cond, __LINE__, __FILE__)
#define compare(val1, val2) compare(val1, val2, __LINE__, __FILE__)
#define compareTolerance(val1, val2, tol) compareTolerance(val1, val2, tol, __LINE__, __FILE__)

using namespace Eigen;

void test_fill_in_k() {
    VectorXd v = VectorXd::Zero(10);
    VectorXd k = VectorXd::Zero(v.size());
    for (size_t i=0;i<v.size();i++) v[i] = i/v.size();
    int p = 1; //eys
    double k_min = 0; 
    double k_max = 1;
    fill_in_k(k,v,k_max,k_min,p);
    tf::compare(true,v.isApprox(k));
    

    tf::check(true);
    tf::check(1); // 1 is equal to true in c++
    tf::compare(1, 1); // Arguments should be of the exact same type
    tf::compareTolerance(1., 1.1, 0.2);

    Eigen::Vector3d v1;
    v1 << 1, 2, 3;
    Eigen::Vector3d v2;
    v2 << 1, 2, 3;

    tf::compare(true, v1.isApprox(v2));
    tf::compareTolerance(0., (v1 - v2).norm(), 0.1);
}

int main() {
    test_fill_in_k();
    return 0;
}