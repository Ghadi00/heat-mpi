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
void run_strips(const Params& p, long steps, double dt, double r)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, P = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);

    const int N = p.N, W = N + 1;
    if (P > N - 1) die("too many processes: strips need N-1 >= P (every rank owns at least one row)");

    const int nloc = split_count(rank, N, P);
    const int j0 = split_first(rank, N, P);   // global row of local row 1
    const int jg = j0 - 1;                    // global row of local row 0 (ghost)
    const int up   = (rank + 1 < P) ? rank + 1 : MPI_PROC_NULL;
    const int down = (rank > 0)     ? rank - 1 : MPI_PROC_NULL;

    // nloc owned rows + 2 ghost rows, N+1 values each.  Fill every row,
    // ghosts included, from the initial condition: ghost rows then hold the
    // neighbour's initial values (refreshed by the exchange anyway), and for
    // rank 0 / P-1 the plate's boundary row, which is never overwritten
    // because the exchange with MPI_PROC_NULL receives nothing.
    const std::size_t nvals = std::size_t(nloc + 2) * W;
    std::vector<double> u_old(nvals), u_new(nvals);
    for (int l = 0; l <= nloc + 1; ++l)
        for (int i = 0; i <= N; ++i)
            u_old[std::size_t(l) * W + i] = u_new[std::size_t(l) * W + i] =
                initial_value(p.test_case, i, jg + l, N);

    auto row = [&](std::vector<double>& u, int l) { return u.data() + std::size_t(l) * W; };

    // Halo exchange on the CURRENT time level u.
    //  shift (default): two shifts in which every rank takes part at once,
    //    1. send my top owned row up,     receive my bottom ghost from below
    //    2. send my bottom owned row down, receive my top ghost from above
    //    so each Sendrecv costs about one message time, whatever P is.
    //  pair: the order written in the project brief,
    //    1. send top row up,     receive top ghost from above (partner: up)
    //    2. send bottom row down, receive bottom ghost from below (partner: down)
    //    Correct, but rank k's first call can only finish once rank k+1 has
    //    reached its second call, so the exchanges ripple down from the top
    //    rank one after another: the cost grows with P.
    auto exchange = [&](std::vector<double>& u) {
        if (p.halo == "shift") {
            MPI_Sendrecv(row(u, nloc), W, MPI_DOUBLE, up, 0,
                         row(u, 0), W, MPI_DOUBLE, down, 0, comm, MPI_STATUS_IGNORE);
            MPI_Sendrecv(row(u, 1), W, MPI_DOUBLE, down, 1,
                         row(u, nloc + 1), W, MPI_DOUBLE, up, 1, comm, MPI_STATUS_IGNORE);
        } else {
            MPI_Sendrecv(row(u, nloc), W, MPI_DOUBLE, up, 2,
                         row(u, nloc + 1), W, MPI_DOUBLE, up, 2, comm, MPI_STATUS_IGNORE);
            MPI_Sendrecv(row(u, 1), W, MPI_DOUBLE, down, 2,
                         row(u, 0), W, MPI_DOUBLE, down, 2, comm, MPI_STATUS_IGNORE);
        }
    };

    Timing mine;
    MPI_Barrier(comm);                 // start the clocks together
    const double t0 = MPI_Wtime();
    for (long n = 0; n < steps; ++n) {
        const double ta = MPI_Wtime();
        exchange(u_old);
        const double tb = MPI_Wtime();
        update_rows(u_old.data(), u_new.data(), nloc, W, r);
        const double tc = MPI_Wtime();
        std::swap(u_old, u_new);
        mine.comm += tb - ta;
        mine.comp += tc - tb;
    }
    mine.total = MPI_Wtime() - t0;

    Timing t;   // maximum over ranks
    MPI_Reduce(&mine.total, &t.total, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&mine.comp,  &t.comp,  1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&mine.comm,  &t.comm,  1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (p.print) {
        exchange(u_old);   // refresh ghosts so they show the current values
        print_in_rank_order(format_rows(u_old, W, nloc + 2, jg, rank, 0, nloc + 1, "ghost"), comm);
    }

    // Errors for case V over the same nodes as the serial version: every
    // owned row, plus the boundary rows j = 0 (rank 0) and j = N (rank P-1).
    double e_inf = 0, e_sum2 = 0;
    const bool have_err = (p.test_case == "V");
    if (have_err) {
        const int jlo = (rank == 0) ? 0 : j0;
        const int jhi = (rank == P - 1) ? N : j0 + nloc - 1;
        double li = 0, ls = 0;
        error_sums(u_old, W, N, 0, 0, N, jg, jlo, jhi, steps * dt, li, ls);
        MPI_Reduce(&li, &e_inf, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&ls, &e_sum2, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
    const double h = 1.0 / N;
    if (rank == 0) print_summary(p, P, steps, r, dt, t, have_err, e_inf, std::sqrt(h * h * e_sum2));

    // Output: gather the owned rows on rank 0.  Rows 0 and N are fixed
    // boundary rows, so rank 0 fills them from the boundary conditions.
    if (!p.out.empty()) {
        std::vector<int> counts, displs;
        std::vector<double> full;
        if (rank == 0) {
            counts.resize(P);
            displs.resize(P);
            for (int k = 0; k < P; ++k) {
                counts[k] = split_count(k, N, P) * W;
                displs[k] = split_first(k, N, P) * W;
            }
            full.assign(std::size_t(W) * W, 0.0);
            for (int i = 0; i <= N; ++i) {
                full[i] = initial_value(p.test_case, i, 0, N);
                full[std::size_t(N) * W + i] = initial_value(p.test_case, i, N, N);
            }
        }
        MPI_Gatherv(row(u_old, 1), nloc * W, MPI_DOUBLE,
                    rank == 0 ? full.data() : nullptr, counts.data(), displs.data(),
                    MPI_DOUBLE, 0, comm);
        if (rank == 0) write_field(p.out, full);
    }
}
}
