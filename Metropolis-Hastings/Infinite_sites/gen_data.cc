#include <cstdlib>
#include <gsl/gsl_spmatrix.h>
#include <iostream>
#include "tree.hh"

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Call " << argv[0] << " <sample size> <mutation rate>" << std::endl;
        return 1;
    }
    int n = atoi(argv[1]);
    double theta = atof(argv[2]);
    Sim sim(n, theta);
    sim.generate_data();
    for (unsigned int i = 0; i < sim.data->size1; i++) {
        for (unsigned int j = 0; j < sim.data->size2; j++) {
            std::cout << gsl_matrix_get(sim.data, i, j) << " ";
        }
        std::cout << sim.row_counts[i] << std::endl;
    }
    return 1;
}
