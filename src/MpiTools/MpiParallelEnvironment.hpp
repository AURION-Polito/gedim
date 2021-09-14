#ifndef __MPIPARALLELENVIRONMENT_H
#define __MPIPARALLELENVIRONMENT_H


#include <memory>
#include "IMpiProcess.hpp"

#if USE_PETSC ==1
#include "petsc.h"
#endif

namespace GeDiM
{
  /// \brief Interface used to implement the Mpi Parallel Environment
  /// \copyright See top level LICENSE file for details.
  class MpiParallelEnvironment
  {
    protected:
      static shared_ptr<IMpiProcess> _process;

    public:
      static const IMpiProcess& Process() { return *_process; }

      /// Initialize the MPI Environment
      /// \param argc command line argument counter
      /// \param argv command live argument values
      template<typename ProcessType>
      static Output::ExitCodes Initialize(int argc,
                                          char** argv)
      {
#if USE_PETSC == 1
        PetscInitialize(&argc, &argv, NULL, NULL);
#elif USE_MPI == 1
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

        return Output::Success;
      }

      /// Finalize the MPI Environment
      static Output::ExitCodes Finalize()
      {
#if USE_PETSC == 1
        PetscFinalize();
#elif USE_MPI == 1
        MPI::Finalize();
#else
        return Output::Success;
#endif
        return Output::Success;
      }
  };
}

#endif // __MPIPARALLELENVIRONMENT_H
