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

#include "MpiProcess.hpp"

#include "IOUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
MpiProcess::MpiProcess(unsigned int rank, unsigned int numberProcesses, bool isActive)
{
    _rank = rank;
    _numberProcesses = numberProcesses;
    _isActive = isActive;
}
MpiProcess::~MpiProcess()
{
}
// ***************************************************************************
void MpiProcess::Summary()
{
    std::string activeLabel = (_isActive ? "Active" : "Not Active");

    Output::PrintLine('-', true);
    Output::PrintGenericMessage("Process %d / %d - %s", false, _rank, _numberProcesses, activeLabel.c_str());
    Output::PrintLine('-', true);
}
// ***************************************************************************
int MpiExtensions::mpiSendTag = 1;
int MpiExtensions::mpiRecvTag = 1;
// ***************************************************************************

} // namespace Gedim
