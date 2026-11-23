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
// Storage (notes, Section 8.1): one contiguous std::vector<double> of
// (N+1)^2 values, row-major: idx(i,j) = j*(N+1) + i, so row j is the line
// y = y_j and i varies fastest.  Two arrays for the two time levels; their
// contents are swapped (std::swap of vectors swaps pointers, no copy).
//
// Output file (--out): the full (N+1)x(N+1) grid, boundary included, raw
// little-endian doubles in row-major order, no header.
//
// Standard output: one summary line of key=value pairs.
//
// Only --decomp serial exists so far; strips (Part 2) and blocks (Part 3)
// will be added to this file.
#include <mpi.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
    bool allow_unstable = false;
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
        "            [--decomp serial] [--allow-unstable]\n");
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

// Grid: (N+1) x (N+1) nodes, row-major.
struct Grid {
    int N;
    int W;   // row width = N + 1
    explicit Grid(int n) : N(n), W(n + 1) {}
    inline std::size_t idx(int i, int j) const { return std::size_t(j) * W + i; }
    std::size_t size() const { return std::size_t(W) * W; }
};

// Initial condition including boundary values, for every node.
double initial_value(const std::string& c, int i, int j, int N)
{
    const bool left = (i == 0), right = (i == N), bottom = (j == 0), top = (j == N);
    const bool boundary = left || right || bottom || top;
    if (c == "1") return boundary ? 1.0 : 0.0;
    if (c == "2") {
        if (bottom || top) return 1.0;     // includes the four corners
        if (left || right) return 0.0;
        return 0.0;
    }
    // "V"
    if (boundary) return 0.0;
    const double h = 1.0 / N;
    return std::sin(kPi * i * h) * std::sin(kPi * j * h);
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

int run_serial(const Params& p, long steps, double dt, double r_used)
{
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) die("--decomp serial must run on exactly one process");

    const Grid g(p.N);
    const int N = p.N;
    std::vector<double> u_old(g.size()), u_new(g.size());

    // Initial condition and boundary values in BOTH arrays: boundary nodes
    // are never updated, so they must already be right in u_new.
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i)
            u_old[g.idx(i, j)] = u_new[g.idx(i, j)] = initial_value(p.test_case, i, j, N);

    const double r = r_used;
    const int W = g.W;

    auto t0 = std::chrono::steady_clock::now();
    for (long n = 0; n < steps; ++n) {
        for (int j = 1; j < N; ++j) {
            const double* uo = &u_old[g.idx(0, j)];
            double* un = &u_new[g.idx(0, j)];
            for (int i = 1; i < N; ++i) {
                // Same expression, same order of operations as the parallel
                // versions will use -> bit-for-bit identical results.
                un[i] = uo[i] + r * (uo[i + 1] + uo[i - 1] + uo[i + W] + uo[i - W] - 4.0 * uo[i]);
            }
        }
        std::swap(u_old, u_new);   // swaps the internal pointers, no copy
    }
    auto t1 = std::chrono::steady_clock::now();
    const double t_total = std::chrono::duration<double>(t1 - t0).count();

    const double T_final = steps * dt;
    std::printf("decomp=serial P=1 case=%s N=%d steps=%ld r=%.6g dt=%.6e T=%.6g "
                "t_total=%.6f t_comp=%.6f t_comm=%.6f",
                p.test_case.c_str(), N, steps, r_used, dt, T_final, t_total, t_total, 0.0);

    if (p.test_case == "V") {
        const double h = 1.0 / N;
        const double decay = std::exp(-2.0 * kPi * kPi * T_final);
        double e_inf = 0.0, e_sum2 = 0.0;
        for (int j = 0; j <= N; ++j)
            for (int i = 0; i <= N; ++i) {
                const double exact = std::sin(kPi * i * h) * std::sin(kPi * j * h) * decay;
                const double e = u_old[g.idx(i, j)] - exact;
                e_inf = std::fmax(e_inf, std::fabs(e));
                e_sum2 += e * e;
            }
        std::printf(" err_inf=%.6e err_2h=%.6e", e_inf, std::sqrt(h * h * e_sum2));
    }
    std::printf("\n");

    if (!p.out.empty()) write_field(p.out, u_old);
    return 0;
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

    int rc = 0;
    if (p.decomp == "serial") rc = run_serial(p, steps, dt, r_used);
    else if (p.decomp == "strips" || p.decomp == "blocks")
        die("--decomp " + p.decomp + " is not implemented yet (Parts 2 and 3)");
    else die("--decomp must be serial, strips or blocks");

    MPI_Finalize();
    return rc;
}
