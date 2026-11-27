#include "verify/error_norms.hpp"
#include "core/constants.hpp"
#include <cmath>
#include <cstddef>

namespace heat {
void error_sums(const std::vector<double>& u, int W, int N,
                int i_of_local0, int ilo, int ihi, int j_of_local0, int jlo, int jhi,
                double T, double& e_inf, double& e_sum2)
{
    const double h = 1.0 / N;
    const double decay = std::exp(-2.0 * kPi * kPi * T);
    for (int j = jlo; j <= jhi; ++j) {
        const double* row = &u[std::size_t(j - j_of_local0) * W];
        for (int i = ilo; i <= ihi; ++i) {
            const double exact = std::sin(kPi * i * h) * std::sin(kPi * j * h) * decay;
            const double e = row[i - i_of_local0] - exact;
            e_inf = std::fmax(e_inf, std::fabs(e));
            e_sum2 += e * e;
        }
    }
}
}
