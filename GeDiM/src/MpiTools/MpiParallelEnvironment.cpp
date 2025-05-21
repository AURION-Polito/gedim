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

#include "MpiParallelEnvironment.hpp"
#include "MpiProcess.hpp"

namespace Gedim
{
// ***************************************************************************
std::shared_ptr<IMpiProcess> MpiParallelEnvironment::_process(new MpiProcess(0, 1, true));
// ***************************************************************************
} // namespace Gedim
