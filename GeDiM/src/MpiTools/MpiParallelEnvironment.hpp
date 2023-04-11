#ifndef __GEDIM_MPIPARALLELENVIRONMENT_H
#define __GEDIM_MPIPARALLELENVIRONMENT_H

#include <memory>
#include "IMpiProcess.hpp"

namespace Gedim
{
  /// \brief Interface used to implement the Mpi Parallel Environment
  /// \copyright See top level LICENSE file for details.
  class MpiParallelEnvironment
  {
    protected:
      static std::shared_ptr<IMpiProcess> _process;

    public:
      static const IMpiProcess& Process() { return *_process; }

      /// Initialize the MPI Environment
      /// \param argc command line argument counter
      /// \param argv command live argument values
      template<typename ProcessType>
      static void Initialize(int argc,
                             char** argv)
      {
#if USE_MPI == 1
        MPI::Init(argc, argv);
#endif

        unsigned int rank = 0;
        unsigned int numberProcesses = 1;

#if USE_MPI == 1
        MPI_Comm_rank(MPI_COMM_WORLD, (int*)&rank);
        MPI_Comm_size(MPI_COMM_WORLD, (int*)&numberProcesses);
#endif // USE_MPI

        _process.reset(new ProcessType(rank,
                                       numberProcesses,
                                       true));
      }

      /// Finalize the MPI Environment
      static void Finalize()
      {
#if USE_MPI == 1
        MPI::Finalize();
#endif
      }
  };
}

#endif // __GEDIM_MPIPARALLELENVIRONMENT_H
