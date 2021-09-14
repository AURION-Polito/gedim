#ifndef __MPIPROCESS_H
#define __MPIPROCESS_H

#include <vector>
#include <string>

#include "Output.hpp"
#include "Eigen"
using namespace Eigen;

#include "IMpiProcess.hpp"

using namespace std;

namespace Gedim
{
  class MpiProcess : public IMpiProcess
  {
    protected:
      unsigned int _rank; ///< Rank of the process
      unsigned int _numberProcesses; ///< Number of all the processes of the MPI simulation
      bool _isActive; ///< Tells if the process is active

    public:
      MpiProcess(unsigned int rank,
                 unsigned int numberProcesses,
                 bool isActive);
      virtual ~MpiProcess();

      unsigned int Rank() const { return _rank; }
      unsigned int NumberProcesses() const { return _numberProcesses; }
      bool IsActive() const { return _isActive; }

      virtual void Summary();
  };

  class MpiExtensions
  {
    private:
      static int mpiSendTag; ///< Tag used for MPI send communications
      static int mpiRecvTag; ///< Tag used for MPI receive communications
    public:
      static int MpiSendTag(const int& numSends = 1) { int tag = mpiSendTag; mpiSendTag += numSends; return tag; }
      static int MpiRecvTag(const int& numRecvs = 1) { int tag = mpiRecvTag; mpiRecvTag += numRecvs; return tag; }
  };
}


#endif // __MPIPROCESS_H

