#ifndef meshSolver_hpp
#define meshSolver_hpp

#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;


Eigen::SparseMatrix<double> find_K(
    const Eigen::VectorXd& v,
    std::vector<std::vector<int>> rectangles,
    int number_of_points,
    const Eigen::Matrix4d& local_matrix,
    double k_min,
    double k_max,
    double penal
) {
    using T = Eigen::Triplet<double>;
    std::vector<T> triplets;

    Eigen::VectorXd k_values = k_min + (k_max - k_min) * v.array().pow(penal);


    int size = number_of_points * number_of_points;
    Eigen::SparseMatrix<double> K(size, size);

    for (size_t e = 0; e < rectangles.size(); ++e) {
        double k_e = k_values[e];

        for (int l = 0; l < 4; ++l) {
            for (int m = 0; m < 4; ++m) {
                int index_i = rectangles[e][l];
                int index_j = rectangles[e][m];
                double value = local_matrix(l, m);

                triplets.emplace_back(index_i, index_j, k_e * value);
            }
        }
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}

Eigen::VectorXd find_F(std::vector<std::vector<int>> rectangles, int number_of_points, double L) {
    double h = (L / 2) / (number_of_points - 1);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(number_of_points * number_of_points);

    for (const auto& rectangle : rectangles) {
        double F_local = (h * h) * (10e6) / 2.;

        for (int l = 0; l < 4; ++l) {
            int index_i = rectangle[l];
            F(index_i) += F_local;
        }
    }

    return F;
}


std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> apply_boundary(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, std::vector<Eigen::Vector3d>& boundary_points, double T_k) {
    int count = 0 ;


    for (Eigen::Vector3d point : boundary_points) {
        int index_of_point = point[0];
        F -= T_k * K.col(index_of_point);  

        for (int i = 0; i < K.rows(); ++i) {
            K.coeffRef(i, index_of_point) = 0;
            K.coeffRef(index_of_point, i) = 0;
        }
        K.coeffRef(index_of_point, index_of_point) = 1; 
        F(index_of_point) = T_k;
        count ++;
    }

    return { K, F };
}


std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> apply_boundary_conditions_band(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<Eigen::Vector3d>& boundary_points, double T_k) {
    int n = K.rows();
    int bandwidth = static_cast<int>(std::sqrt(n)) + 1;  

    for (const auto& point : boundary_points) {
        int index = static_cast<int>(point[0]); 
        F -= T_k * K.col(index); 
    
        for (int k = std::max(0, index - 4*bandwidth); k < std::min(n, index + 4*bandwidth); ++k) {
            if (k != index) {  
                K.coeffRef(k, index) = 0.0;
                K.coeffRef(index, k) = 0.0;
            }
        }
        K.coeffRef(index, index) = 1.0;

        F(index) = T_k;
    }
    return { K, F };
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> apply_boundary_conditions_optimized(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<Eigen::Vector3d>& boundary_points, double T_k) {
    int n = K.rows();
    int bandwidth = static_cast<int>(std::sqrt(n)) + 1;  

    for (const auto& point : boundary_points) {
        int index = static_cast<int>(point[0]); 

        K.coeffRef(index, index) = 1e31 ;

        F(index) = T_k*1e31;

    }
    return { K, F };
}



#endif