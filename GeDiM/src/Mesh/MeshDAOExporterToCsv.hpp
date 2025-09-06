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

#ifndef __MeshDAOExporterToCsv_H
#define __MeshDAOExporterToCsv_H

#include "IMeshDAO.hpp"
#include "MeshFromCsvUtilities.hpp"

namespace Gedim
{
/// \brief MeshDAOExporterToCsv
/// \copyright See top level LICENSE file for details.
class MeshDAOExporterToCsv final
{
  private:
    const MeshFromCsvUtilities &utilities;

  public:
    MeshDAOExporterToCsv(const MeshFromCsvUtilities &utilities);
    ~MeshDAOExporterToCsv();

    /// \brief Export the mesh in all parts
    /// \param configuration the configuration for export
    /// \param mesh the mesh to be exported
    void Export(const MeshFromCsvUtilities::Configuration &configuration, const IMeshDAO &mesh) const;
};

} // namespace Gedim

#endif // __MeshDAOExporterToCsv_H
