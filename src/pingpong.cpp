// pingpong.cpp -- measure point-to-point latency and bandwidth between two
// MPI processes (Part 0, Step 6 of the project).
//
// Rank 0 sends m bytes to rank 1, rank 1 sends them straight back. This is
// repeated many times for every message size m = 8, 16, ..., max bytes.
// Rank 0 prints one CSV line per size:
//
//   bytes, reps, rtt_avg_s, rtt_median_s, rtt_min_s, oneway_avg_s, bandwidth_MBps
//
// The one-way time is half the average round trip; it is what the model
// t_msg(m) = t_s + m * t_w is fitted to (see scripts/fit_pingpong.py).
//
// The number of repetitions adapts to the network: for each size we time two
// warm-up round trips and choose reps so that the size takes about --budget
// seconds (clamped to [--min-reps, --max-reps]). A fixed 1000 reps would take
// minutes per size over a 200 ms Tailscale link.
//
// Options:
//   --max BYTES      largest message size (default 16777216 = 16 MiB)
//   --budget SEC     target measuring time per size (default 2.0)
//   --min-reps N     default 5
//   --max-reps N     default 1000
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank = 0, size = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (size != 2) {
        if (rank == 0) std::fprintf(stderr, "pingpong needs exactly 2 processes (got %d)\n", size);
        MPI_Abort(comm, 1);
    }

    long max_bytes = 16L * 1024 * 1024;
    double budget = 2.0;
    int min_reps = 5, max_reps = 1000;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) { if (rank == 0) std::fprintf(stderr, "missing value for %s\n", a.c_str()); MPI_Abort(comm, 1); }
            return argv[++i];
        };
        if      (a == "--max")      max_bytes = std::atol(next());
        else if (a == "--budget")   budget    = std::atof(next());
        else if (a == "--min-reps") min_reps  = std::atoi(next());
        else if (a == "--max-reps") max_reps  = std::atoi(next());
        else { if (rank == 0) std::fprintf(stderr, "unknown option %s\n", a.c_str()); MPI_Abort(comm, 1); }
    }

    const int peer = 1 - rank;

    // Exchange host names so the output records which pair was measured.
    char me[MPI_MAX_PROCESSOR_NAME] = {0}, other[MPI_MAX_PROCESSOR_NAME] = {0};
    int len = 0;
    MPI_Get_processor_name(me, &len);
    MPI_Sendrecv(me, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, peer, 0,
                 other, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, peer, 0,
                 comm, MPI_STATUS_IGNORE);

    std::vector<char> buf(static_cast<size_t>(max_bytes), 1);

    if (rank == 0) {
        std::printf("# pingpong %s <-> %s\n", me, other);
        std::printf("bytes,reps,rtt_avg_s,rtt_median_s,rtt_min_s,oneway_avg_s,bandwidth_MBps\n");
        std::fflush(stdout);
    }

    for (long m = 8; m <= max_bytes; m *= 2) {
        const int count = static_cast<int>(m);

        auto round_trip = [&]() {
            if (rank == 0) {
                MPI_Send(buf.data(), count, MPI_BYTE, peer, 1, comm);
                MPI_Recv(buf.data(), count, MPI_BYTE, peer, 1, comm, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(buf.data(), count, MPI_BYTE, peer, 1, comm, MPI_STATUS_IGNORE);
                MPI_Send(buf.data(), count, MPI_BYTE, peer, 1, comm);
            }
        };

        // Warm-up (also opens the connection for the first size) and estimate.
        double t0 = MPI_Wtime();
        round_trip();
        round_trip();
        double estimate = (MPI_Wtime() - t0) / 2.0;

        int reps = 0;
        if (rank == 0)
            reps = std::clamp(static_cast<int>(budget / std::max(estimate, 1e-9)), min_reps, max_reps);
        MPI_Bcast(&reps, 1, MPI_INT, 0, comm);

        std::vector<double> rtt(reps);
        MPI_Barrier(comm);
        for (int r = 0; r < reps; ++r) {
            double t = MPI_Wtime();
            round_trip();
            rtt[r] = MPI_Wtime() - t;
        }

        if (rank == 0) {
            double sum = 0.0;
            for (double x : rtt) sum += x;
            double avg = sum / reps;
            std::vector<double> sorted = rtt;
            std::sort(sorted.begin(), sorted.end());
            double median = (reps % 2) ? sorted[reps / 2]
                                       : 0.5 * (sorted[reps / 2 - 1] + sorted[reps / 2]);
            double oneway = avg / 2.0;
            std::printf("%ld,%d,%.9e,%.9e,%.9e,%.9e,%.3f\n",
                        m, reps, avg, median, sorted.front(), oneway, m / oneway / 1e6);
            std::fflush(stdout);
        }
    }

    MPI_Finalize();
    return 0;
}
