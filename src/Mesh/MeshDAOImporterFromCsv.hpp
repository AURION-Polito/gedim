#ifndef __MeshDAOImporterFromCsv_H
#define __MeshDAOImporterFromCsv_H

#include "IMeshDAO.hpp"
#include "MeshImporterFromCsvUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  /// \brief MeshDAOImporterFromCsv
  /// \note each file could be EmptyFileReader if not necessary
  /// \copyright See top level LICENSE file for details
  class MeshDAOImporterFromCsv final {
    public:
      MeshDAOImporterFromCsv();
      ~MeshDAOImporterFromCsv();

      void Import(MeshImporterFromCsvUtilities::Configuration& configuration,
                  MeshImporterFromCsvUtilities& importerUtilities,
                  IMeshDAO& mesh);

      void ImportMesh2D(MeshImporterFromCsvUtilities::Configuration& configuration,
                        MeshImporterFromCsvUtilities& importerUtilities,
                        IMeshDAO& mesh);
  };

}

#endif // __MeshDAOImporterFromCsv_H
