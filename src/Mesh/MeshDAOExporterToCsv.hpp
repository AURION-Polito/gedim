#ifndef __MeshDAOExporterToCsv_H
#define __MeshDAOExporterToCsv_H

#include "IMeshDAO.hpp"
#include "MeshFromCsvUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  /// \brief MeshDAOExporterToCsv
  /// \copyright See top level LICENSE file for details.
  class MeshDAOExporterToCsv final
  {
    private:
      const MeshFromCsvUtilities& utilities;

    public:
      MeshDAOExporterToCsv(const MeshFromCsvUtilities& utilities);
      ~MeshDAOExporterToCsv();

      /// \brief Export the mesh in all parts
      /// \param configuration the configuration for export
      /// \param mesh the mesh to be exported
      void Export(const MeshFromCsvUtilities::Configuration& configuration,
                  const IMeshDAO& mesh) const;
  };

}

#endif // __MeshDAOExporterToCsv_H
