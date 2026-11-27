#pragma once
#include <algorithm>
namespace heat {
// Even split of the N-1 interior nodes of one direction over `parts` parts:
// part k owns nodes first(k) .. first(k)+count(k)-1; the first (N-1) mod parts get one extra.
inline int split_count(int k, int N, int parts) { return (N - 1) / parts + (k < (N - 1) % parts ? 1 : 0); }
inline int split_first(int k, int N, int parts) { return 1 + k * ((N - 1) / parts) + std::min(k, (N - 1) % parts); }
}
