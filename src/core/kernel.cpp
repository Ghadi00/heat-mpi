#include "core/kernel.hpp"
#include <cstddef>

namespace heat {
void update_rows(const double* uo, double* un, int nrows, int W, double r)
{
    for (int j = 1; j <= nrows; ++j) {
        const double* a = uo + std::size_t(j) * W;
        double* b = un + std::size_t(j) * W;
        for (int i = 1; i < W - 1; ++i)
            b[i] = a[i] + r * (a[i + 1] + a[i - 1] + a[i + W] + a[i - W] - 4.0 * a[i]);
    }
}
}
