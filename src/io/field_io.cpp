#include "io/field_io.hpp"
#include "utilities/errors.hpp"
#include <cstdint>
#include <cstdio>

namespace heat {
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
}
