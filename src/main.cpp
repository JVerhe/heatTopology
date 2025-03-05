#include <Eigen/Sparse>
#include <iostream>
#include "optimization.hpp"
#include "meshHelper.hpp"


void readConfig(const std::string& filename, int& number_of_points, int& p, int& ft) {
    std::ifstream file(filename);
    file >> number_of_points >> p >> ft;
    file.close();
}

void writeOutput(const std::string& filename, Eigen::VectorXd x){
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i = 0; i < x.size(); ++i) {
            file << x(i) << std::endl;
        }
        file.close();
        std::cout << "Saved to " << filename << std::endl;
    }
    else {
        std::cerr << "Error when saving " << filename << std::endl;
    }
    return;

}


int main() {

    int number_of_points; int p; int ft;
    char config_file[100] = "../Config_files/config.txt";

    readConfig(config_file,number_of_points,p,ft);

    double  L = 0.01;
    double T_k = 293;
    double vol_frac = 0.4;

    Eigen::Matrix4d local_matrix;
    local_matrix << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;


    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(number_of_points);
    Eigen::VectorXd x = optimize(local_matrix, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, L, T_k, ft);


    char outputfile[100] = "../Results/results.txt";
    save_result_to_file(x, outputfile);

    return 0;
}