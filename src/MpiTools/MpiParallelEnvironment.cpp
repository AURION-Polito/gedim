#include "MpiParallelEnvironment.hpp"
#include "MpiProcess.hpp"

namespace Gedim
{
  // ***************************************************************************
  shared_ptr<IMpiProcess> MpiParallelEnvironment::_process(new MpiProcess(0, 1, true));
  // ***************************************************************************
}
