#include "mpi.h"

#include <array>
#include <cassert>
#include <cstring> // for memcpy
#include <iostream>
#include <iterator>
#include <numeric>
#include <ostream>
#include <vector>

extern "C"
{
//  void __sll_m_collective_MOD_sll_s_boot_collective(int32_t *mpi_mode);
  void sll_s_allocate_collective();
  void sll_s_set_communicator_collective(MPI_Fint *mpi_comm);
  void sll_s_halt_collective();

  void sim_bsl_vp_3d3v_cart_dd_slim_init(void** sim, const char *filename);
  void sim_bsl_vp_3d3v_cart_dd_slim_run(void** sim);
  void sim_bsl_vp_3d3v_cart_dd_slim_delete(void** sim);
  void sim_bsl_vp_3d3v_cart_dd_slim_get_distribution(void** sim, double *& cPtr);
  void sim_bsl_vp_3d3v_cart_dd_slim_set_distribution(void** sim, double *& cPtr);
  void sim_bsl_vp_3d3v_cart_dd_slim_get_local_size(void** sim, int32_t *cPtr);
  void sim_bsl_vp_3d3v_cart_dd_slim_advect_v(void** sim, double *delta_t);
  void sim_bsl_vp_3d3v_cart_dd_slim_advect_x(void** sim, double *delta_t);
  void sim_bsl_vp_3d3v_cart_dd_slim_print_etas(void** sim);
  void sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init(void** sim);
  void sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics(void** sim, int32_t* timeStepNumber);
}

// helper function to output any array
template<typename T, size_t N>
inline std::ostream &operator <<(std::ostream &os, const std::array<T, N> &v) {
  using namespace std;
  os << "[";
  copy(v.begin(), v.end(), ostream_iterator<T>(os, ", "));
  os << "]";
  return os;
}

// helper function to output any vector
template<typename T>
inline std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
  using namespace std;
  os << "[";
  copy(v.begin(), v.end(), ostream_iterator<T>(os, ", "));
  os << "]";
  return os;
}

void initMpiSelalibStyle(int argc, char **argv, int mpi_mode = MPI_THREAD_FUNNELED)
{
  sll_s_allocate_collective();
  int ignore;
  MPI_Init_thread(&argc, &argv, mpi_mode, &ignore);
}

int main(int argc, char **argv)
{
#ifdef _OPENMP
  int mpi_mode = MPI_THREAD_MULTIPLE;
#else
  int mpi_mode = MPI_THREAD_SINGLE;
#endif

  MPI_Comm mpi_world = MPI_COMM_WORLD;

  // this would call MPI_Init in selalib:
  // __sll_m_collective_MOD_sll_s_boot_collective(&mpi_mode);
  // replaced by the following:
  initMpiSelalibStyle(argc, argv, mpi_mode);

  auto f_mpi_world = MPI_Comm_c2f(mpi_world);
  int rank = 0;
  MPI_Comm_rank(mpi_world, &rank);
  sll_s_set_communicator_collective(&f_mpi_world);

  // to get two separate communicators
  // int colour = rank % 2;
  // to get the same communicator
  int colour = 0;
  MPI_Comm newComm;
  MPI_Comm_split(mpi_world, colour, 0, &newComm);
  auto f_newComm = MPI_Comm_c2f(newComm);
  sll_s_set_communicator_collective(&f_newComm);

  void* simPtr;
  auto simPtrPtr = &simPtr;
  std::cout << "init" << std::endl;
  assert(argc > 1);
  sim_bsl_vp_3d3v_cart_dd_slim_init(simPtrPtr, argv[1]);

  sim_bsl_vp_3d3v_cart_dd_slim_print_etas(simPtrPtr);

  std::cout << "get local size" << std::endl;
  std::array<int32_t, 6> localSize;
  sim_bsl_vp_3d3v_cart_dd_slim_get_local_size(simPtrPtr, localSize.data());
  std::cout << localSize << std::endl;

  bool doDiagnostics = true;
  if (doDiagnostics) {
    // have a "diagnostics-only" simulation that does only diagnostics output (which should stay the same for all time steps)
    sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init(simPtrPtr);

    for (int i=1; i<4; ++i){
      int32_t* iPtr = &i;
      sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics(simPtrPtr,iPtr);
    }
  } else{

    int32_t bufferSize = std::accumulate(localSize.begin(), localSize.end(), 1, std::multiplies<int32_t>());

    std::vector<double> field (bufferSize, 1.);
    for (int i=0; i<3; ++i){
      std::cout << "get distribution" << std::endl;
      double * fieldAddress;
      sim_bsl_vp_3d3v_cart_dd_slim_get_distribution(simPtrPtr, fieldAddress);
      field = std::vector<double>(fieldAddress, fieldAddress+bufferSize);
      //std::cout << field << std::endl;
      // std::cout << std::accumulate(field.begin(), field.end(), 0) << std::endl;
      // // modify field
      // for (auto & f: field){
      //   f += 1.;
      // }
      // std::cout << std::accumulate(field.begin(), field.end(), 0) << std::endl;
      auto fieldData = field.data();
      std::cout << "set distribution" << std::endl;
      // sim_bsl_vp_3d3v_cart_dd_slim_set_distribution(sim, fieldData);
      std::memcpy(fieldAddress, fieldData, bufferSize*sizeof(double));

      std::cout << "run" << std::endl;
      sim_bsl_vp_3d3v_cart_dd_slim_run(simPtrPtr);
    }
  }

  std::cout << "delete" << std::endl;
  sim_bsl_vp_3d3v_cart_dd_slim_delete(simPtrPtr);

  std::cout << "works in cpp" << std::endl;

  sll_s_halt_collective();
}
