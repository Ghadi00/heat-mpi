#pragma once
#include <string>
namespace heat {
// Initial condition including boundary values, for global node (i, j).
double initial_value(const std::string& c, int i, int j, int N);
}
