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
    cfg.readFile(argv[1]);
    std::string input_file;
    cfg.lookupValue("input_file", input_file);
    Sim sim(input_file);
    int run_length;
    cfg.lookupValue("run_length", run_length);
    double sd_mutation_rate;
    double sd_merger_times;
    cfg.lookupValue("sd_mutation_rate", sd_mutation_rate);
    cfg.lookupValue("sd_merger_times", sd_merger_times);
    int accept_topology = 0;
    int accept_times = 0;
    int accept_mutation = 0;
    for (int i = 0; i < run_length; i++) {
        accept_topology += sim.metropolis_topology();
        accept_times += sim.metropolis_times(sd_merger_times);
        accept_mutation += sim.metropolis_mutation_rate(sd_mutation_rate);
        std::cout << 2 * sim.mutation_rate << " " << sim.tree_height() << " "
                  << std::endl;
    }
    std::cerr << "acceptance probabilities:" << std::endl;
    std::cerr << "mutation rate:\t" << (double)accept_mutation / run_length
              << std::endl;
    std::cerr << "merger times:\t" << (double)accept_times / run_length
              << std::endl;
    std::cerr << "subtree-prune-regraft:\t"
              << (double)accept_topology / run_length << std::endl;
    return 0;
}
