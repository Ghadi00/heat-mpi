#pragma once
#include "cli/params.hpp"
#include "common/timing.hpp"
#include <mpi.h>
#include <string>
#include <vector>
namespace heat {
// Local rows (top first) with their global row index, for --print.
std::string format_rows(const std::vector<double>& u, int W, int nrows_total, int j_of_local0,
                        int rank, int ghost_lo, int ghost_hi, const char* edge_tag,
                        int i_of_local0 = 0, const std::string& title = "");
// Every rank's text printed by rank 0, in rank order.
void print_in_rank_order(const std::string& mine, MPI_Comm comm);
// One key=value summary line (call on rank 0 only).
void print_summary(const Params& p, int P, long steps, double r_used, double dt,
                   const Timing& t, bool have_err, double e_inf, double e_2h,
                   const std::string& extra = "");
}
