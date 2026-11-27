#pragma once
#include <string>
namespace heat {
struct Params {
    std::string test_case = "1";   // "1", "2" or "V"
    int N = 64;
    double r = 0.2;                // requested alpha*dt/h^2
    double T = -1.0;               // final time (if > 0)
    long steps = -1;               // fixed number of steps (if >= 0)
    std::string out;               // output file, empty = none
    std::string decomp = "serial";
    std::string halo = "shift";    // strips: "shift" (default) or "pair"
    bool allow_unstable = false;
    bool print = false;            // debug: print local arrays incl. ghosts
};
void usage();
Params parse(int argc, char** argv);
}
