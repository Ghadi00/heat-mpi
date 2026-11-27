#pragma once
namespace heat {
// THE update rule (Task 5.2), shared by every decomposition.
// Rows 1..nrows of `un` are computed from `uo`; rows 0 and nrows+1 are only read.
void update_rows(const double* uo, double* un, int nrows, int W, double r);
}
