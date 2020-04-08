#include <cstdlib>
#include <iostream>
#include "tree.hh"

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Call " << argv[0] << " <path to config file>" << std::endl;
        return 1;
    }
    Sim sim(argv[1]);
    double accept_spr = 0;
    double accept_theta = 0;
    sim.run(accept_spr, accept_theta);
    std::cout << accept_spr << " " << accept_theta << std::endl;
    return 1;
}
