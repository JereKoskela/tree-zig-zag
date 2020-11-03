#include "tree.hh"
#include <cstdlib>
#include <iostream>
#include <libconfig.h++>

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Call " << argv[0] << " <config file>" << std::endl;
        return 1;
    }
    libconfig::Config cfg;
    std::string input_file;
    double run_length, mutation_velocity, hybrid_rate, hybrid_mutation_rate_sd;
    double spr_acceptance_rate = 0;
    double mutation_rate_acceptance_rate = 0;

    cfg.readFile(argv[1]);
    cfg.lookupValue("input_file", input_file);
    cfg.lookupValue("run_length", run_length);
    cfg.lookupValue("mutation_velocity", mutation_velocity);
    cfg.lookupValue("hybrid_rate", hybrid_rate);
    cfg.lookupValue("hybrid_mutation_rate_sd", hybrid_mutation_rate_sd);

    Sim sim(input_file, mutation_velocity, hybrid_rate,
            hybrid_mutation_rate_sd);
    std::cout << sim.mutation_rate << " " << sim.tree_height();
    if (hybrid_rate > 0) {
        std::cout << " 0";
    }
    std::cout << std::endl;
    sim.run(run_length, spr_acceptance_rate, mutation_rate_acceptance_rate);
    if (hybrid_rate > 0) {
        std::cerr << "Metropolis acceptance probabilities:" << std::endl;
        std::cerr << "mutation rate:\t" << mutation_rate_acceptance_rate
                  << std::endl;
        std::cerr << "subtree-prune-regraft:\t" << spr_acceptance_rate
                  << std::endl;
    }
    return 0;
}
