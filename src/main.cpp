#include <Eigen/Sparse>
#include <iostream>
#include <fstream>


/**
 * Creates a list of rectangles, where each rectangle is defined by the indices of its four nodes.
 * The output is a vector of size (number_of_rectangles × 4), where each row represents a rectangle.
 * 
 * @param number_of_points The number of mesh points in one fourth of the mesh.
 * @return A vector containing lists of four node indices forming each rectangle.
 */
std::vector<std::vector<int>> create_rectangle_and_mesh(int number_of_points) {

    Eigen::MatrixXi mesh = Eigen::MatrixXi::Zero(number_of_points, number_of_points);
    
    int count = 0;
    for (int i = 0; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points; ++j) {
            mesh(i, j) = count++;
        }
    }

    std::vector<std::vector<int>> rectangles;
    for (int i = 1; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points - 1; ++j) {
            std::vector<int> line = {
                mesh(i, j),
                mesh(i, j + 1),
                mesh(i - 1, j + 1),
                mesh(i - 1, j)
            };
            rectangles.push_back(line);
        }
    }

    return rectangles;
}


/**
 * Generates a list of coordinates for mesh points in one-fourth of the domain.
 * The output is a matrix of size (number_of_points² × 2), where each row represents a (x, y) coordinate.
 *
 * @param L The length of the chip.
 * @param number_of_points The number of mesh points in one fourth of the mesh.
 * @return An Eigen matrix of size (number_of_points² × 2) containing the coordinates.
 */
Eigen::MatrixXd create_coordinates(double L, int number_of_points) {
    Eigen::MatrixXd coordinates(number_of_points * number_of_points, 2);
    double h = (L / 2.0) / (number_of_points - 1);

    for (int i = 0; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points; ++j) {
            coordinates(i * number_of_points + j, 0) = j * h; // x-coordinate
            coordinates(i * number_of_points + j, 1) = i * h; // y-coordinate
        }
    }

    return coordinates;
}

std::vector<Eigen::Vector3d> filter_boundary_points_with_index(const Eigen::MatrixXd& coordinates, double L, double lower = 0.003, double upper = 0.007) {
    /**
     * Creates a list of boundary points with their indices.
     * Each boundary point is stored as a 3D vector (index, x, y).
     * 
     * @param coordinates: Matrix (N x 2) containing the (x, y) coordinates.
     * @param L: The length of the chip.
     * @param lower: Lower bound for y values.
     * @param upper: Upper bound for y values.
     * @return A vector of Eigen::Vector3d where each entry contains (index, x, y).
     */
    
    std::vector<Eigen::Vector3d> boundary_points;
    
    for (int idx = 0; idx < coordinates.rows(); ++idx) {
        double x = coordinates(idx, 0);
        double y = coordinates(idx, 1);

        if ((x == 0 || x == L) && (lower <= y && y <= upper)) {
            boundary_points.emplace_back(idx, x, y);
        }
    }

    return boundary_points;
}

Eigen::VectorXd fill_in_k(const Eigen::VectorXd& v, double k_max, double k_min, double p) {
    /**
     * Computes the material property values (k) based on the density vector (v).
     * 
     * @param v: Vector of densities (values between 0 and 1).
     * @param k_max: Maximum stiffness value.
     * @param k_min: Minimum stiffness value.
     * @param p: Penalization factor.
     * @return Eigen::VectorXd containing computed k values.
     */

    Eigen::VectorXd k(v.size());

    for (int i = 0; i < v.size(); ++i) {
        k(i) = k_min + (k_max - k_min) * std::pow(v(i), p);
    }

    return k;
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

    return {K, F};
}


void save_result_to_file(const Eigen::VectorXd& U, const std::string& filename) {
    std::ofstream file(filename);  // Ouvre un fichier pour écrire
    if (file.is_open()) {
        for (int i = 0; i < U.size(); ++i) {
            file << U(i) << std::endl;  // Écrit chaque valeur de U dans le fichier
        }
        file.close();  // Ferme le fichier après l'écriture
        std::cout << "Saved" << filename << std::endl;
    } else {
        std::cerr << "Not saved"<< filename << std::endl;
    }
}




int main() {
    double  L = 0.01;
    int p = 3;
    double T_k  = 293;
    int number_of_points = 70;

    Eigen::Matrix4d local_matrix; // Matrice 4x4 de type double

    // Remplissage de la matrice
    local_matrix <<  2.0/3, -1.0/6, -1.0/3, -1.0/6,
                    -1.0/6,  2.0/3, -1.0/6, -1.0/3,
                    -1.0/3, -1.0/6,  2.0/3, -1.0/6,
                    -1.0/6, -1.0/3, -1.0/6,  2.0/3;




    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(number_of_points);
    Eigen::VectorXd v = Eigen::VectorXd::Ones(rectangles.size()) * 0.2;
    Eigen::MatrixXd coordinates = create_coordinates(L,number_of_points);
    std::vector<Eigen::Vector3d> boundary_points = filter_boundary_points_with_index(coordinates,L);
    Eigen::SparseMatrix<double> K = find_K(v, rectangles, number_of_points, local_matrix, 0.2, 65.0, p);
    Eigen::VectorXd F = find_F(rectangles,number_of_points,L);
    auto result = apply_boundary(K, F, boundary_points, T_k);
    K = result.first;
    F = result.second;


    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
    solver.compute(K);  // Initialisation du solver

    // Application d'un préconditionneur, si disponible
    solver.setMaxIterations(1000);
    solver.setTolerance(1e-6);

    Eigen::VectorXd U = solver.solve(F);  // Résolution du système

    if (solver.error() > solver.tolerance()) {
        std::cerr << "La résolution a échoué!" << std::endl;
        return -1;
    }

    std::cout << "Solution U : " << std::endl;

    save_result_to_file(U,"proto/result.txt");

    
    return 0;
}