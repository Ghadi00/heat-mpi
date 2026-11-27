#include "utilities/errors.hpp"
#include <mpi.h>
#include <cstdio>
#include <cstdlib>

namespace heat {
void die(const std::string& msg)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::fprintf(stderr, "heat: %s\n", msg.c_str());
    MPI_Finalize();
    std::exit(1);
}

void fatal(const std::string& msg)
{
    std::fprintf(stderr, "heat: %s\n", msg.c_str());
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::exit(1);
}
}
