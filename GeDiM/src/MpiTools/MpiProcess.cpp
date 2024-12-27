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
