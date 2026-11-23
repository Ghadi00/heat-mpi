// heat.cpp -- explicit finite-difference solver for the 2D heat equation
//
//     u_t = alpha (u_xx + u_yy)  on [0,1]^2,  alpha = 1,
//
// on a uniform grid of N x N intervals (h = 1/N, nodes 0 <= i,j <= N).
//
// Scheme (background notes, Task 5.2, with dx = dy = h and r = alpha dt / h^2):
//
//   u_new(i,j) = u(i,j) + r * ( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4 u(i,j) )
//
// Stability (Task 6.3): r <= 1/4.  The program refuses to run above that
// unless --allow-unstable is given (needed only for Task 6.7).
//
// Cases (project brief, Section 3):
//   1  all four edges held at u = 1, interior starts at 0
//   2  top (y=1) and bottom (y=0) edges at u = 1 incl. corners, left/right at 0
//   V  verification: u0 = sin(pi x) sin(pi y), u = 0 on all edges,
//      exact solution sin(pi x) sin(pi y) exp(-2 pi^2 t)
//
// Storage (notes, Section 8.1): contiguous std::vector<double>, row-major,
// row j is the line y = y_j and i varies fastest; a row has W = N+1 values.
// Two arrays for the two time levels, swapped (std::swap of vectors swaps
// pointers, no copy).
//
// Decompositions (--decomp):
//   serial  one process, the whole (N+1)x(N+1) grid.  The reference.
//   strips  (Part 2) the N-1 interior rows are split into P contiguous strips;
//           the first (N-1) mod P ranks get one extra row; rank 0 owns the
//           bottom strip.  Each rank stores its nloc rows plus one ghost row
//           below (local row 0) and above (local row nloc+1); for rank 0 and
//           P-1 the ghost row is the plate's boundary row.  Ghost rows are
//           refreshed by a halo exchange with MPI_Sendrecv before every step.
//   blocks  (Part 3) not implemented yet.
//
// Serial and strips call the SAME update function (update_rows), so every
// node is computed by identical arithmetic in identical order and the output
// files are bit-for-bit identical for any number of processes.
//
// Output file (--out): the full (N+1)x(N+1) grid, boundary included, raw
// little-endian doubles, row-major, no header.  In parallel, rank 0 gathers
// the rows with MPI_Gatherv and writes the file.  Without --out nothing is
// gathered (for large N over a slow network that would take minutes).
//
// Standard output: one summary line of key=value pairs, printed by rank 0.
// t_total, t_comp, t_comm are each the maximum over ranks of the time spent
// in the time-stepping loop, in its computation, and in its halo exchange.
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRMax = 0.25;   // stability limit from Task 6.3 (2D, dx = dy)

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

// Error that every rank detects identically (bad command line): print once,
// shut MPI down cleanly, exit with status 1.
[[noreturn]] void die(const std::string& msg)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::fprintf(stderr, "heat: %s\n", msg.c_str());
    MPI_Finalize();
    std::exit(1);
}

// Error on one rank only (e.g. cannot write a file): abort all ranks.
[[noreturn]] void fatal(const std::string& msg)
{
    std::fprintf(stderr, "heat: %s\n", msg.c_str());
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::exit(1);
}

void usage()
{
    std::fprintf(stderr,
        "usage: heat --case 1|2|V --N <intervals> --r <alpha dt/h^2>\n"
        "            (--T <final time> | --steps <n>) [--out file.bin]\n"
        "            [--decomp serial|strips] [--halo shift|pair]\n"
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

// Initial condition including boundary values, for global node (i, j).
double initial_value(const std::string& c, int i, int j, int N)
{
    const bool left = (i == 0), right = (i == N), bottom = (j == 0), top = (j == N);
    const bool boundary = left || right || bottom || top;
    if (c == "1") return boundary ? 1.0 : 0.0;
    if (c == "2") {
        if (bottom || top) return 1.0;     // includes the four corners
        return 0.0;                        // left/right edges and interior
    }
    // "V"
    if (boundary) return 0.0;
    const double h = 1.0 / N;
    return std::sin(kPi * i * h) * std::sin(kPi * j * h);
}

// THE update rule (Task 5.2), shared by every decomposition.
// `uo` and `un` point at local row 0 of an array whose rows have W values.
// Rows 1..nrows are updated; rows 0 and nrows+1 (ghost or boundary rows) are
// only read.  Columns 0 and W-1 are the left/right boundary: not updated.
void update_rows(const double* uo, double* un, int nrows, int W, double r)
{
    for (int j = 1; j <= nrows; ++j) {
        const double* a = uo + std::size_t(j) * W;
        double* b = un + std::size_t(j) * W;
        for (int i = 1; i < W - 1; ++i)
            b[i] = a[i] + r * (a[i + 1] + a[i - 1] + a[i + W] + a[i - W] - 4.0 * a[i]);
    }
}

void write_field(const std::string& path, const std::vector<double>& u)
{
    // The file format is little-endian; x86-64 and ARM64 Linux both are.
    const uint16_t probe = 1;
    if (*reinterpret_cast<const unsigned char*>(&probe) != 1)
        fatal("big-endian machine: byte swapping not implemented");
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) fatal("cannot open " + path + " for writing");
    if (std::fwrite(u.data(), sizeof(double), u.size(), f) != u.size())
        fatal("short write to " + path);
    std::fclose(f);
}

// Error against the exact solution of case V, accumulated over global rows
// jlo..jhi.  Local row l of `u` is global row j_of_local0 + l.
void error_sums(const std::vector<double>& u, int W, int j_of_local0, int jlo, int jhi,
                double T, double& e_inf, double& e_sum2)
{
    const int N = W - 1;
    const double h = 1.0 / N;
    const double decay = std::exp(-2.0 * kPi * kPi * T);
    for (int j = jlo; j <= jhi; ++j) {
        const double* row = &u[std::size_t(j - j_of_local0) * W];
        for (int i = 0; i <= N; ++i) {
            const double exact = std::sin(kPi * i * h) * std::sin(kPi * j * h) * decay;
            const double e = row[i] - exact;
            e_inf = std::fmax(e_inf, std::fabs(e));
            e_sum2 += e * e;
        }
    }
}

// Format local rows (top first, like the plate) with their global row index.
std::string format_rows(const std::vector<double>& u, int W, int nrows_total, int j_of_local0,
                        int rank, int ghost_lo, int ghost_hi, const char* edge_tag)
{
    std::string s;
    char buf[128];
    std::snprintf(buf, sizeof buf, "---- rank %d: local rows 0..%d = global rows %d..%d\n",
                  rank, nrows_total - 1, j_of_local0, j_of_local0 + nrows_total - 1);
    s += buf;
    for (int l = nrows_total - 1; l >= 0; --l) {
        const char* tag = (l == ghost_lo || l == ghost_hi) ? edge_tag : "     ";
        std::snprintf(buf, sizeof buf, "  %s l=%2d j=%2d |", tag, l, j_of_local0 + l);
        s += buf;
        for (int i = 0; i < W; ++i) {
            std::snprintf(buf, sizeof buf, " %7.4f", u[std::size_t(l) * W + i]);
            s += buf;
        }
        s += "\n";
    }
    return s;
}

// Print every rank's text in rank order.  Output from different processes
// can arrive interleaved, so the texts are sent to rank 0, which prints them.
void print_in_rank_order(const std::string& mine, MPI_Comm comm)
{
    int rank = 0, P = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);
    int len = static_cast<int>(mine.size());
    std::vector<int> lens(P), displs(P);
    MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, 0, comm);
    std::string all;
    if (rank == 0) {
        int total = 0;
        for (int k = 0; k < P; ++k) { displs[k] = total; total += lens[k]; }
        all.resize(total);
    }
    MPI_Gatherv(mine.data(), len, MPI_CHAR, rank == 0 ? &all[0] : nullptr,
                lens.data(), displs.data(), MPI_CHAR, 0, comm);
    if (rank == 0) { std::fputs(all.c_str(), stdout); std::fflush(stdout); }
}

struct Timing { double total = 0, comp = 0, comm = 0; };

void print_summary(const Params& p, int P, long steps, double r_used, double dt,
                   const Timing& t, bool have_err, double e_inf, double e_2h)
{
    std::printf("decomp=%s P=%d N=%d steps=%ld t_total=%.6f t_comp=%.6f t_comm=%.6f "
                "case=%s r=%.6g dt=%.6e T=%.6g",
                p.decomp.c_str(), P, p.N, steps, t.total, t.comp, t.comm,
                p.test_case.c_str(), r_used, dt, steps * dt);
    if (p.decomp == "strips") std::printf(" halo=%s", p.halo.c_str());
    if (have_err) std::printf(" err_inf=%.6e err_2h=%.6e", e_inf, e_2h);
    std::printf("\n");
    std::fflush(stdout);
}

// ---------------------------------------------------------------- serial --
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
    if (have_err) error_sums(u_old, W, 0, 0, N, steps * dt, e_inf, e_sum2);
    const double h = 1.0 / N;
    print_summary(p, 1, steps, r, dt, t, have_err, e_inf, std::sqrt(h * h * e_sum2));

    if (!p.out.empty()) write_field(p.out, u_old);
}

// ---------------------------------------------------------------- strips --
// Rows owned by rank k: global rows first(k) .. first(k)+count(k)-1.
int strip_count(int k, int N, int P) { return (N - 1) / P + (k < (N - 1) % P ? 1 : 0); }
int strip_first(int k, int N, int P) { return 1 + k * ((N - 1) / P) + std::min(k, (N - 1) % P); }

void run_strips(const Params& p, long steps, double dt, double r)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, P = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);

    const int N = p.N, W = N + 1;
    if (P > N - 1) die("too many processes: strips need N-1 >= P (every rank owns at least one row)");

    const int nloc = strip_count(rank, N, P);
    const int j0 = strip_first(rank, N, P);   // global row of local row 1
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
        error_sums(u_old, W, jg, jlo, jhi, steps * dt, li, ls);
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
                counts[k] = strip_count(k, N, P) * W;
                displs[k] = strip_first(k, N, P) * W;
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

}  // namespace

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
    else if (p.decomp == "blocks") die("--decomp blocks is not implemented yet (Part 3)");
    else die("--decomp must be serial, strips or blocks");

    MPI_Finalize();
    return 0;
}
