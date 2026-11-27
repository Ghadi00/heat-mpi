#include "cli/params.hpp"
#include "common/errors.hpp"
#include "core/constants.hpp"
#include <mpi.h>
#include <cstdio>
#include <cstdlib>

namespace heat {
void usage()
{
    std::fprintf(stderr,
        "usage: heat --case 1|2|V --N <intervals> --r <alpha dt/h^2>\n"
        "            (--T <final time> | --steps <n>) [--out file.bin]\n"
        "            [--decomp serial|strips|blocks] [--halo shift|pair]\n"
        "            [--print] [--allow-unstable]\n");
}

Params parse(int argc, char** argv)
{
    Params p;
    for (int k = 1; k < argc; ++k) {
        std::string a = argv[k];
        auto val = [&]() -> std::string {
            if (k + 1 >= argc) die("missing value for " + a);
            return argv[++k];
        };
        if      (a == "--case")   p.test_case = val();
        else if (a == "--N")      p.N = std::atoi(val().c_str());
        else if (a == "--r")      p.r = std::atof(val().c_str());
        else if (a == "--T")      p.T = std::atof(val().c_str());
        else if (a == "--steps")  p.steps = std::atol(val().c_str());
        else if (a == "--out")    p.out = val();
        else if (a == "--decomp") p.decomp = val();
        else if (a == "--halo")   p.halo = val();
        else if (a == "--print")  p.print = true;
        else if (a == "--allow-unstable") p.allow_unstable = true;
        else if (a == "--help" || a == "-h") { usage(); MPI_Finalize(); std::exit(0); }
        else { usage(); die("unknown option " + a); }
    }
    if (p.test_case == "v") p.test_case = "V";
    if (p.test_case != "1" && p.test_case != "2" && p.test_case != "V")
        die("--case must be 1, 2 or V");
    if (p.N < 2) die("--N must be at least 2");
    if (!(p.r > 0)) die("--r must be positive");
    if ((p.T > 0) == (p.steps >= 0)) die("give exactly one of --T or --steps");
    if (p.halo != "shift" && p.halo != "pair") die("--halo must be shift or pair");
    if (p.r > kRMax && !p.allow_unstable) {
        char buf[256];
        std::snprintf(buf, sizeof buf,
            "r = %g is unstable: the explicit 2D scheme needs r = alpha*dt/h^2 <= 1/4 "
            "(von Neumann analysis, Task 6.3). Use a smaller --r "
            "(or --allow-unstable to watch it blow up, Task 6.7).", p.r);
        die(buf);
    }
    return p;
}
}
