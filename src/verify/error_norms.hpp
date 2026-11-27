#pragma once
#include <vector>
namespace heat {
// Error against the exact solution of case V over global nodes ilo..ihi x jlo..jhi.
// Local node (l, m) of `u` (row length W) is global node (i_of_local0 + l, j_of_local0 + m).
void error_sums(const std::vector<double>& u, int W, int N,
                int i_of_local0, int ilo, int ihi, int j_of_local0, int jlo, int jhi,
                double T, double& e_inf, double& e_sum2);
}
