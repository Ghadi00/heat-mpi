#include "solvers/solvers.hpp"
#include "utilities/errors.hpp"
#include "utilities/timing.hpp"
#include "core/initial.hpp"
#include "core/kernel.hpp"
#include "core/partition.hpp"
#include "io/field_io.hpp"
#include "io/report.hpp"
#include "verify/error_norms.hpp"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace heat {
void run_serial(const Params& p, long steps, double dt, double r)
{
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) die("--decomp serial must run on exactly one process");

    const int N = p.N, W = N + 1;
    std::vector<double> u_old(std::size_t(W) * W), u_new(std::size_t(W) * W);

    // Initial condition and boundary values in BOTH arrays: boundary nodes
    // are never updated, so they must already be right in u_new.
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i)
            u_old[std::size_t(j) * W + i] = u_new[std::size_t(j) * W + i] =
                initial_value(p.test_case, i, j, N);

    // The whole grid is one "strip" of N-1 interior rows whose ghost rows
    // are the bottom and top boundary rows.
    const double t0 = MPI_Wtime();
    for (long n = 0; n < steps; ++n) {
        update_rows(u_old.data(), u_new.data(), N - 1, W, r);
        std::swap(u_old, u_new);
    }
    Timing t;
    t.total = t.comp = MPI_Wtime() - t0;

    if (p.print) print_in_rank_order(format_rows(u_old, W, N + 1, 0, 0, 0, N, "bound"), MPI_COMM_WORLD);

    double e_inf = 0, e_sum2 = 0;
    const bool have_err = (p.test_case == "V");
    if (have_err) error_sums(u_old, W, N, 0, 0, N, 0, 0, N, steps * dt, e_inf, e_sum2);
    const double h = 1.0 / N;
    print_summary(p, 1, steps, r, dt, t, have_err, e_inf, std::sqrt(h * h * e_sum2));

    if (!p.out.empty()) write_field(p.out, u_old);
}
}
