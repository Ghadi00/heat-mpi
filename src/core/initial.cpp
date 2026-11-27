#include "core/initial.hpp"
#include "core/constants.hpp"
#include <cmath>

namespace heat {
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
}
