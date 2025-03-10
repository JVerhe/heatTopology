#include <iostream>
#include <Eigen/Sparse>
#include <cassert>

// Include the next lines in your own files
#include <testFramework.hpp>
#include <../src/optHelper.hpp>
#include <../src/meshHelper.hpp>
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
    tf::compareTolerance(0.,(v-k).norm(),0.001);

    p=2;
    fill_in_k(k,v,k_max,k_min,p);
    VectorXd k_compare = VectorXd::Zero(v.size());
    for (size_t i=0;i<v.size();i++) k_compare[i] = (i/v.size())*(i/v.size());
    tf::compareTolerance(0.,(k_compare-k).norm(),0.001);
    

    tf::check(true);
    tf::check(1); // 1 is equal to true in c++
    tf::compare(1, 1); // Arguments should be of the exact same type
    tf::compareTolerance(1., 1.1, 0.2);

    Eigen::Vector3d v1;
    v1 << 1, 2, 3;
    Eigen::Vector3d v2;
    v2 << 1, 2, 3;
    tf::compareTolerance(0., (v1 - v2).norm(), 0.1);
}

void test_objective() {
    int p = 1; 
    double k_min = 0; 
    double k_max = 1;
    double N_elements_1d = 3;
    int N_points_1d = N_elements_1d+1;
    double T_const = 300;
    VectorXd v = VectorXd::Ones(N_elements_1d*N_elements_1d);
    Eigen::Matrix4d K0;
    K0 << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;
    const Eigen::VectorXd T = T_const*VectorXd::Ones(N_points_1d*N_points_1d);
    
    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(N_points_1d);

    double obj = objective(v, rectangles, T, K0, k_min,k_max,p);
    const Eigen::VectorXd T_e = T_const*VectorXd::Ones(4);
    Eigen::VectorXd K0T = (N_elements_1d*N_elements_1d) * (K0*T_e);
    double obj_ref = T_e.transpose() * K0T;
    
    tf::compareTolerance(obj, obj_ref, 0.0001);
}

void test_update_densities() {
    int N_total_elements = 4;
    double move = 0.1;
    VectorXd x_old = VectorXd::Ones(N_total_elements);
    VectorXd x_new = VectorXd::Zero(N_total_elements);
    VectorXd B = VectorXd::Zero(N_total_elements);
    B[0] = 0; B[1] = 1.1; B[2] = 0.95; B[3] = 0.92;

    update_densities(x_old, x_new, B, move);
    Eigen::VectorXd x_new_ref = VectorXd::Zero(N_total_elements);
    x_new_ref[0] = 0.9; x_new_ref[1]=1;x_new_ref[2] = 0.95;x_new_ref[3] = 0.92;

    tf::compareTolerance((x_new-x_new_ref).norm(),0.,0.0001);
}



int main() {
    test_fill_in_k();
    test_objective();
    test_update_densities();
    return 0;
}