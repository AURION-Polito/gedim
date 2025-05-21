// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __GEDIM_IMPIPROCESS_H
#define __GEDIM_IMPIPROCESS_H

namespace Gedim
{
/// \brief Interface used to implement the MPI process
/// \copyright See top level LICENSE file for details.
class IMpiProcess
{
  public:
    virtual ~IMpiProcess()
    {
    }

    /// \return the rank of the MPI process
    virtual unsigned int Rank() const = 0;
    /// \return the total number of processes in the communicator
    virtual unsigned int NumberProcesses() const = 0;
    /// \return if the MPI process is active
    virtual bool IsActive() const = 0;

    /// \brief Summary of the process
    virtual void Summary() = 0;
};
} // namespace Gedim

#endif // __GEDIM_IMPIPROCESS_H
