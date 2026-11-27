#pragma once
#include "cli/params.hpp"
namespace heat {
void run_serial(const Params& p, long steps, double dt, double r);
void run_strips(const Params& p, long steps, double dt, double r);
void run_blocks(const Params& p, long steps, double dt, double r);
}
