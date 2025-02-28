#ifndef meshSolver_hpp
#define meshSolver_hpp

#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;

void createSparseMatrix(int nx, int ny, const vector<int>& iH, const vector<int>& jH, const vector<double>& sH, SparseMatrix<double>& H, VectorXd& Hs) {
    int N = nx * ny;
    vector<Triplet<double>> tripletList;

    for (size_t k = 0; k < sH.size(); ++k) {
        tripletList.emplace_back(iH[k], jH[k], sH[k]);
    }

    H.resize(N, N);
    H.setFromTriplets(tripletList.begin(), tripletList.end());

    // Compute sum of each row
    Hs = VectorXd::Zero(N);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            Hs(it.row()) += it.value();
        }
    }
}

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
        double F_local = (h * h) * (10e6) / 4.;

        for (int l = 0; l < 4; ++l) {
            int index_i = rectangle[l];
            F(index_i) += F_local;
        }
    }

    return F;
}


std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> apply_boundary(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, std::vector<Eigen::Vector3d> boundary_points, double T_k) {
    for (const auto& point : boundary_points) {
        int index_of_point = point[0];
        F -= T_k * K.col(index_of_point);  // Soustraction du terme de force local

        // Mise à zéro des coefficients de la matrice K associés au point de bordure
        for (int i = 0; i < K.rows(); ++i) {
            K.coeffRef(i, index_of_point) = 0;
            K.coeffRef(index_of_point, i) = 0;
        }

        K.coeffRef(index_of_point, index_of_point) = 1; // Diagonal égale à 1
        F(index_of_point) = T_k; // Force au point de bordure
    }

    return { K, F };
}
#endif