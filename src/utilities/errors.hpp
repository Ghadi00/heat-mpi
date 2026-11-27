#pragma once
#include <string>

namespace heat {
// Error every rank detects identically (bad command line): print once, finalize MPI, exit(1).
[[noreturn]] void die(const std::string& msg);
// Error on one rank only (e.g. cannot write a file): abort all ranks.
[[noreturn]] void fatal(const std::string& msg);
}
