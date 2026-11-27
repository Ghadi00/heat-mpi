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
void run_blocks(const Params& p, long steps, double dt, double r)
{
    int P = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    const int N = p.N;

    // Process grid: dims[0] = Px (x, columns), dims[1] = Py (y, rows).
    int dims[2] = {0, 0};
    MPI_Dims_create(P, 2, dims);
    const int Px = dims[0], Py = dims[1];
    if (Px > N - 1 || Py > N - 1)
        die("too many processes: blocks need N-1 >= Px and N-1 >= Py");

    int periods[2] = {0, 0};          // non-periodic: the plate has edges
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, /*reorder=*/0, &cart);
    int rank = 0, coords[2] = {0, 0};
    MPI_Comm_rank(cart, &rank);
    MPI_Cart_coords(cart, rank, 2, coords);
    const int cx = coords[0], cy = coords[1];

    int left, right, down, up;        // MPI_PROC_NULL outside the plate
    MPI_Cart_shift(cart, 0, 1, &left, &right);
    MPI_Cart_shift(cart, 1, 1, &down, &up);

    const int nx = split_count(cx, N, Px), i0 = split_first(cx, N, Px);   // my columns
    const int ny = split_count(cy, N, Py), j0 = split_first(cy, N, Py);   // my rows
    const int W = nx + 2;             // local row length incl. left/right ghosts
    const int ig = i0 - 1, jg = j0 - 1;   // global node of local (0, 0)

    // (ny+2) x (nx+2) local array.  Fill everything, ghosts included, from
    // the initial condition: ghosts on the plate's edge then hold the
    // boundary values and are never overwritten (MPI_PROC_NULL receives
    // nothing); the other ghosts are refreshed by every exchange.
    const std::size_t nvals = std::size_t(ny + 2) * W;
    std::vector<double> u_old(nvals), u_new(nvals);
    for (int m = 0; m <= ny + 1; ++m)
        for (int l = 0; l <= nx + 1; ++l)
            u_old[std::size_t(m) * W + l] = u_new[std::size_t(m) * W + l] =
                initial_value(p.test_case, ig + l, jg + m, N);

    // A column of my block: ny values, one per row, W doubles apart.
    MPI_Datatype column;
    MPI_Type_vector(ny, 1, W, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    auto at = [&](std::vector<double>& u, int l, int m) { return u.data() + std::size_t(m) * W + l; };

    // Four shifts, like the strips' default: everyone sends up, then down,
    // then right, then left.  Rows: nx contiguous doubles (columns 1..nx).
    // Columns: one `column`, starting in row 1.  Corners are not sent.
    auto exchange = [&](std::vector<double>& u) {
        MPI_Sendrecv(at(u, 1, ny), nx, MPI_DOUBLE, up, 0,
                     at(u, 1, 0), nx, MPI_DOUBLE, down, 0, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv(at(u, 1, 1), nx, MPI_DOUBLE, down, 1,
                     at(u, 1, ny + 1), nx, MPI_DOUBLE, up, 1, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv(at(u, nx, 1), 1, column, right, 2,
                     at(u, 0, 1), 1, column, left, 2, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv(at(u, 1, 1), 1, column, left, 3,
                     at(u, nx + 1, 1), 1, column, right, 3, cart, MPI_STATUS_IGNORE);
    };

    Timing mine;
    MPI_Barrier(cart);
    const double t0 = MPI_Wtime();
    for (long n = 0; n < steps; ++n) {
        const double ta = MPI_Wtime();
        exchange(u_old);
        const double tb = MPI_Wtime();
        update_rows(u_old.data(), u_new.data(), ny, W, r);   // same kernel as serial
        const double tc = MPI_Wtime();
        std::swap(u_old, u_new);
        mine.comm += tb - ta;
        mine.comp += tc - tb;
    }
    mine.total = MPI_Wtime() - t0;

    Timing t;
    MPI_Reduce(&mine.total, &t.total, 1, MPI_DOUBLE, MPI_MAX, 0, cart);
    MPI_Reduce(&mine.comp,  &t.comp,  1, MPI_DOUBLE, MPI_MAX, 0, cart);
    MPI_Reduce(&mine.comm,  &t.comm,  1, MPI_DOUBLE, MPI_MAX, 0, cart);

    if (p.print) {
        exchange(u_old);
        char title[96];
        std::snprintf(title, sizeof title, " (coords %d,%d of %dx%d), columns %d..%d",
                      cx, cy, Px, Py, ig, ig + nx + 1);
        print_in_rank_order(format_rows(u_old, W, ny + 2, jg, rank, 0, ny + 1, "ghost", ig, title), cart);
    }

    // Case V errors over owned nodes, plus the plate's boundary nodes held
    // in the ghost layers of the blocks on the edge (so the same nodes as
    // the serial version are counted).
    double e_inf = 0, e_sum2 = 0;
    const bool have_err = (p.test_case == "V");
    if (have_err) {
        const int ilo = (cx == 0) ? 0 : i0, ihi = (cx == Px - 1) ? N : i0 + nx - 1;
        const int jlo = (cy == 0) ? 0 : j0, jhi = (cy == Py - 1) ? N : j0 + ny - 1;
        double li = 0, ls = 0;
        error_sums(u_old, W, N, ig, ilo, ihi, jg, jlo, jhi, steps * dt, li, ls);
        MPI_Reduce(&li, &e_inf, 1, MPI_DOUBLE, MPI_MAX, 0, cart);
        MPI_Reduce(&ls, &e_sum2, 1, MPI_DOUBLE, MPI_SUM, 0, cart);
    }
    const double h = 1.0 / N;
    char extra[32];
    std::snprintf(extra, sizeof extra, "Px=%d Py=%d", Px, Py);
    if (rank == 0) print_summary(p, P, steps, r, dt, t, have_err, e_inf, std::sqrt(h * h * e_sum2), extra);

    // Output: a block is not contiguous in the global grid, so each rank
    // copies its owned nodes into a contiguous buffer, rank 0 gathers the
    // buffers and puts each block in place using that rank's coordinates.
    // The plate's boundary nodes never change: rank 0 fills them itself.
    if (!p.out.empty()) {
        std::vector<double> mine_buf(std::size_t(nx) * ny);
        for (int m = 1; m <= ny; ++m)
            for (int l = 1; l <= nx; ++l)
                mine_buf[std::size_t(m - 1) * nx + (l - 1)] = u_old[std::size_t(m) * W + l];

        std::vector<int> counts, displs;
        std::vector<double> all;
        if (rank == 0) {
            counts.resize(P);
            displs.resize(P);
            int total = 0;
            for (int k = 0; k < P; ++k) {
                int c[2];
                MPI_Cart_coords(cart, k, 2, c);
                counts[k] = split_count(c[0], N, Px) * split_count(c[1], N, Py);
                displs[k] = total;
                total += counts[k];
            }
            all.resize(total);
        }
        MPI_Gatherv(mine_buf.data(), nx * ny, MPI_DOUBLE,
                    rank == 0 ? all.data() : nullptr, counts.data(), displs.data(),
                    MPI_DOUBLE, 0, cart);
        if (rank == 0) {
            const int Wg = N + 1;
            std::vector<double> full(std::size_t(Wg) * Wg);
            for (int j = 0; j <= N; ++j)
                for (int i = 0; i <= N; ++i)
                    if (i == 0 || j == 0 || i == N || j == N)
                        full[std::size_t(j) * Wg + i] = initial_value(p.test_case, i, j, N);
            for (int k = 0; k < P; ++k) {
                int c[2];
                MPI_Cart_coords(cart, k, 2, c);
                const int kx = split_count(c[0], N, Px), ki = split_first(c[0], N, Px);
                const int ky = split_count(c[1], N, Py), kj = split_first(c[1], N, Py);
                const double* src = all.data() + displs[k];
                for (int m = 0; m < ky; ++m)
                    for (int l = 0; l < kx; ++l)
                        full[std::size_t(kj + m) * Wg + (ki + l)] = src[std::size_t(m) * kx + l];
            }
            write_field(p.out, full);
        }
    }

    MPI_Type_free(&column);
    MPI_Comm_free(&cart);
}
}
