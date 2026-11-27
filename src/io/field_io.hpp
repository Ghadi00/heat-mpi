#pragma once
#include <string>
#include <vector>
namespace heat {
// Raw little-endian doubles, row-major, no header.
void write_field(const std::string& path, const std::vector<double>& u);
}
