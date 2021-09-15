#ifndef __GEDIM_IMPIPROCESS_H
#define __GEDIM_IMPIPROCESS_H

#include "IOUtilities.hpp"

using namespace std;

namespace Gedim
{
  /// \brief Interface used to implement the MPI process
  /// \copyright See top level LICENSE file for details.
  class IMpiProcess
	{
		public:
      virtual ~IMpiProcess() { }

      /// \return the rank of the MPI process
      virtual unsigned int Rank() const = 0;
      /// \return the total number of processes in the communicator
      virtual unsigned int NumberProcesses() const = 0;
      /// \return if the MPI process is active
      virtual bool IsActive() const = 0;

      /// \brief Summary of the process
      virtual void Summary() = 0;

	};
}

#endif // __GEDIM_IMPIPROCESS_H

