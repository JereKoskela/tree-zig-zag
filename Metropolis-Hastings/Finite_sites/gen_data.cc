#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <iostream>
#include "tree.hh"

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cout << "Call " << argv[0] << " <sample size> <number of sites> <mutation rate>" << std::endl;
        return 1;
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    double theta = atof(argv[3]);
    Sim sim(n, m, theta);
    sim.generate_data();
    for (unsigned int i = 0; i < sim.data->size1; i++) {
        for (int j = 0; j < m; j++) {
            std::cout << gsl_matrix_get(sim.data, i, j) << " ";
        }
        std::cout << sim.row_counts[i] << std::endl;
    }
    return 1;
}
