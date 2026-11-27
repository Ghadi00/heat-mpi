#include "cli/params.hpp"
#include "common/errors.hpp"
#include "solvers/solvers.hpp"
#include <mpi.h>
#include <cmath>

using namespace heat;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Params p = parse(argc, argv);

    // Time step.  With --T, choose the largest dt <= r h^2 such that T/dt is
    // an integer; the r actually used may then be slightly below the request.
    const double h = 1.0 / p.N;
    long steps;
    double dt;
    if (p.T > 0) {
        const double dt_max = p.r * h * h;
        steps = static_cast<long>(std::ceil(p.T / dt_max - 1e-9));
        dt = p.T / steps;
    } else {
        steps = p.steps;
        dt = p.r * h * h;
    }
    const double r_used = dt / (h * h);

    if (p.decomp == "serial")      run_serial(p, steps, dt, r_used);
    else if (p.decomp == "strips") run_strips(p, steps, dt, r_used);
    else if (p.decomp == "blocks") run_blocks(p, steps, dt, r_used);
    else die("--decomp must be serial, strips or blocks");

    MPI_Finalize();
    return 0;
}
