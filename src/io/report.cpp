#include "io/report.hpp"
#include <cstdio>

namespace heat {
std::string format_rows(const std::vector<double>& u, int W, int nrows_total, int j_of_local0,
                        int rank, int ghost_lo, int ghost_hi, const char* edge_tag,
                        int i_of_local0, const std::string& title)
{
    std::string s;
    char buf[160];
    std::snprintf(buf, sizeof buf, "---- rank %d%s: local rows 0..%d = global rows %d..%d\n",
                  rank, title.c_str(), nrows_total - 1, j_of_local0, j_of_local0 + nrows_total - 1);
    s += buf;
    s += "                   i =";
    for (int i = 0; i < W; ++i) {
        std::snprintf(buf, sizeof buf, " %7d", i_of_local0 + i);
        s += buf;
    }
    s += "\n";
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

void print_summary(const Params& p, int P, long steps, double r_used, double dt,
                   const Timing& t, bool have_err, double e_inf, double e_2h,
                   const std::string& extra)
{
    std::printf("decomp=%s P=%d N=%d steps=%ld t_total=%.6f t_comp=%.6f t_comm=%.6f "
                "case=%s r=%.6g dt=%.6e T=%.6g",
                p.decomp.c_str(), P, p.N, steps, t.total, t.comp, t.comm,
                p.test_case.c_str(), r_used, dt, steps * dt);
    if (p.decomp == "strips") std::printf(" halo=%s", p.halo.c_str());
    if (!extra.empty()) std::printf(" %s", extra.c_str());
    if (have_err) std::printf(" err_inf=%.6e err_2h=%.6e", e_inf, e_2h);
    std::printf("\n");
    std::fflush(stdout);
}
}
