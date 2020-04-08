#include <cstdlib>
#include <iostream>
#include "tree.hh"

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Call " << argv[0] << " <path to config file>" << std::endl;
        return 1;
    }
    Sim sim(argv[1]);
    sim.run();
    return 1;
}
