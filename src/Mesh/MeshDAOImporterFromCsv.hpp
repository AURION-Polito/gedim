#ifndef __MeshDAOImporterFromCsv_H
#define __MeshDAOImporterFromCsv_H

#include "IMeshDAO.hpp"
#include "MeshFromCsvUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  /// \brief MeshDAOImporterFromCsv
  /// \note each file could be EmptyFileReader if not necessary
  /// \copyright See top level LICENSE file for details
  class MeshDAOImporterFromCsv final
  {
    private:
      const MeshFromCsvUtilities& utilities;

    public:
      MeshDAOImporterFromCsv(const MeshFromCsvUtilities& utilities);
      ~MeshDAOImporterFromCsv();

      void Import(MeshFromCsvUtilities::Configuration& configuration,
                  IMeshDAO& mesh);

      void ImportMesh2D(MeshFromCsvUtilities::Configuration& configuration,
                        IMeshDAO& mesh);
  };

}

#endif // __MeshDAOImporterFromCsv_H
