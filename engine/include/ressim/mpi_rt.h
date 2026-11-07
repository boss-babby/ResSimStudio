// C16: MPI loaded at run time (no link-time dependency, the app starts without MPI installed).
//   Windows: msmpi.dll (Microsoft MPI); Linux: libmpi.so.12 / libmpich.so.12 (MPICH, Intel MPI, MVAPICH).
// All of these share the MPICH ABI (integer handles, MPI_COMM_WORLD = 0x44000000, 5-int MPI_Status), so one set of
// constants serves them; Open MPI (pointer handles) is detected and rejected with a message.
// Used by ResSimCLI: "mpiexec -n N ResSimCLI --mpi-batch <jobdir>" distributes study / ensemble runs over the ranks
// of one or many machines (rank 0 dispatches cases dynamically, every rank simulates, results return as messages).
#pragma once
#include <string>

namespace ressim {
namespace mpi {

bool load(std::string* error = nullptr);       // find and bind the library (idempotent)
bool loaded();
std::string libraryName();                     // the library that was loaded
bool init();                                   // MPI_Init (after load)
void finalize();
int rank();
int size();
void barrier();
void abort(int code);
std::string processorName();
// byte messages (any size up to 2 GB): tag >= 0; source -1 = any source, tag -1 = any tag
bool send(const std::string& data, int dest, int tag);
// wait for a matching message without burning a core (short spin, then Iprobe + 0.1-1 ms sleeps); true when one is pending
bool politeWait(int source, int tag);
bool recv(std::string& data, int source, int tag, int* srcOut = nullptr, int* tagOut = nullptr);
double allreduceSum(double v);
double allreduceMax(double v);
// C18 collectives / typed point-to-point for the domain-decomposed linear solver
bool allreduceSumN(double* v, int n);                      // in place
bool bcast(void* data, int bytes, int root);
bool sendDoubles(const double* v, int n, int dest, int tag);
bool recvDoubles(double* v, int n, int source, int tag);
bool sendInts(const int* v, int n, int dest, int tag);
bool recvInts(int* v, int n, int source, int tag);
// deadlock-free pairwise exchange (MPI_Sendrecv)
bool sendrecvDoubles(const double* s, int ns, int dest, double* r, int nr, int source, int tag);
bool sendrecvInts(const int* s, int ns, int dest, int* r, int nr, int source, int tag);

// distributed batch over the ranks: jobdir holds case_XXXX.DATA|.rsim (+ optional job.txt); rank 0 writes
// case_XXXX.rsres (summary + statistics; title "OK" or "ERROR: ...") and progress.txt ("done total"), then done.txt
int runBatchJob(const std::string& jobdir, bool verbose);

} // namespace mpi
} // namespace ressim
