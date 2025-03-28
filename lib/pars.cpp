#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

struct InvParams {
    int Ns;
    int Nt;
    int n_points;
    int trash;
    std::string boot_corr_file;
    int Nboot;
    std::string mean_corr_file;
    double Estar;
    double apar;
    double dLim;
    double uLim;
    std::string out_file;
    std::string rho_vs_lambda_file;
    std::string smear_delta_rec_file;
    std::string rs_vs_lambda;
    double sigma;
    double E0;
};

void read_input(const std::string& input_file_name, InvParams& param) {
    std::ifstream input_fp(input_file_name);
    if (!input_fp) {
        std::cerr << "Error opening file: " << input_file_name << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(input_fp, line)) {
        std::istringstream stream(line);
        std::string key;
        stream >> key;

        if (key == "size") {
            stream >> param.Ns >> param.Nt;
        }
        else if (key == "points_for_inversion") {
            stream >> param.n_points;
        }
        else if (key == "trash") {
            stream >> param.trash;
        }
        else if (key == "boostrap_corr_file") {
            stream >> param.boot_corr_file;
        }
        else if (key == "n_boot") {
            stream >> param.Nboot;
        }
        else if (key == "mean_corr_file") {
            stream >> param.mean_corr_file;
        }
        else if (key == "omega") {
            stream >> param.Estar;
        }
        else if (key == "apar") {
            stream >> param.apar;
        }
        else if (key == "Lside") {
            stream >> param.dLim;
        }
        else if (key == "Rside") {
            stream >> param.uLim;
        }
        else if (key == "output_file") {
            stream >> param.out_file;
        }
        else if (key == "rho_vs_lambda_file") {
            stream >> param.rho_vs_lambda_file;
        }
        else if (key == "smear_func_file") {
            stream >> param.smear_delta_rec_file;
        }
        else if (key == "rs_vs_lambda_file") {
            stream >> param.rs_vs_lambda;
        }
        else if (key == "sigma") {
            std::string temp;
            stream >> temp;
            param.sigma = std::stod(temp);
        }
        else if (key == "Ezero") {
            std::string temp;
            stream >> temp;
            param.E0 = std::stod(temp);
        }
        else {
            std::cerr << "Error: unrecognized option " << key << " in file: " << input_file_name << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    input_fp.close();
}
