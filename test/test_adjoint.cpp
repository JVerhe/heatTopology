#include <iostream>
#include <Eigen/Sparse>
#include <cassert>

// Include the next lines in your own files
#include <testFramework.hpp>
#include <../src/meshHelper.hpp>
#include <../src/optHelper.hpp>
#include <../src/optimization.hpp>
#include <../src/meshSolver.hpp>
#define check(cond) check(cond, __LINE__, __FILE__)
#define compare(val1, val2) compare(val1, val2, __LINE__, __FILE__)
#define compareTolerance(val1, val2, tol) compareTolerance(val1, val2, tol, __LINE__, __FILE__)

using namespace Eigen;

void test_adjoint() 
{
    Eigen::Matrix4d K0;
    K0 << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;
    int nx = 5;
    int p = 2;
    double L = 0.01;

    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(nx+1);
    VectorXd grad = VectorXd::Zero(nx*nx);
    VectorXd v = VectorXd::Constant(nx*nx,0.4);

    MatrixXd coordinates = create_coordinates(L, nx+1);
    std::vector<Vector3d> boundary_points = filter_boundary_points_with_index(coordinates, L);
    SparseMatrix<double> K = find_K(v, rectangles, nx+1, K0, 0.2, 65.0, p);
    VectorXd F = find_F(rectangles, nx+1, L);
    auto result = apply_boundary(K, F, boundary_points, 293);
    K = result.first;
    F = result.second;
    VectorXd U = VectorXd::Zero(F.size());
    solve_sparse_lin_sys(K, F, U);

    const double old_obj = objective(v, rectangles, U, K0, 0.2, 65, p);
    adjoint(grad,U,v,rectangles,K0,p);

    //std::cout << "grad=" << grad << std::endl;
    VectorXd eps_vec = VectorXd::Zero(32);
    VectorXd mean_error = VectorXd::Zero(32);
    VectorXd mean_square_error = VectorXd::Zero(32);
    VectorXd max_error = VectorXd::Zero(32);

    for (int j=0;j<32;j++){
        double eps = std::pow(10,-(j/2)-1);
        if (j%2 != 0) eps*=0.5;
        eps_vec[j] = eps;
    
        VectorXd Fin_diff = VectorXd::Zero(grad.size());
        VectorXd Expected_diff = VectorXd::Zero(grad.size());
        
        //adjust v a tiny bit and compute new objective
        for (size_t i=0;i<grad.size();i++){
            VectorXd v_new = v;
            v_new[i]+=eps;

            SparseMatrix<double> K_new = find_K(v_new, rectangles, nx+1, K0, 0.2, 65.0, p);
            VectorXd F_new = find_F(rectangles, nx+1, L);
            auto result = apply_boundary(K_new, F_new, boundary_points, 293);
            K_new = result.first;
            F_new = result.second;
            VectorXd U_new = VectorXd::Zero(F_new.size());
            solve_sparse_lin_sys(K_new, F_new, U_new);
            double new_obj = objective(v_new, rectangles, U_new, K0, 0.2, 65.0, p);

            Fin_diff[i] = (new_obj-old_obj)*eps; //expect it to be <0
            Expected_diff[i] = grad[i]*eps; //should be <0
            assert(Expected_diff[i]<=0);
        }
        //std::cout << "difference : " << Fin_diff << std::endl;
        //std::cout << "expected difference : " << Expected_diff << std::endl;
        //std::cout << "error : " << (Fin_diff - Expected_diff).cwiseAbs() << std::endl;
        mean_error[j] = (Fin_diff - Expected_diff).cwiseAbs().mean();
        mean_square_error[j] = mean_error[j]*mean_error[j];
        max_error[j] = (Fin_diff - Expected_diff).cwiseAbs().maxCoeff();
    }
    std::cout << "eps_vec : " << eps_vec << std::endl;
    std::cout << "mean_errors : " << mean_error << std::endl;
    std::cout << "max_errors : " << max_error << std::endl;

    save_result_to_file(eps_vec,"adjoint_eps_vec.txt");
    save_result_to_file(mean_error,"adjoint_mean_error.txt");
    save_result_to_file(mean_square_error,"adjoint_MSE.txt");
}

int main() {
    test_adjoint();
    return 0;
}