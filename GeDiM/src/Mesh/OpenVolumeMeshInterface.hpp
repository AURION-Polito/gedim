#ifndef __OpenVolumeMeshInterface_H
#define __OpenVolumeMeshInterface_H

#include "IMeshDAO.hpp"

namespace Gedim
{
  class OpenVolumeMeshInterface final
  {
    public:
      struct OVMMesh final
      {
          struct Cell3D final
          {
              std::vector<unsigned int> FacesIndex;
              std::vector<bool> FacesOrientation;
          };

          unsigned int NumCell0Ds;
          unsigned int NumCell1Ds;
          unsigned int NumCell2Ds;
          unsigned int NumCell3Ds;

          Eigen::MatrixXd Cell0Ds;
          Eigen::MatrixXi Cell1Ds;
          std::vector<Eigen::MatrixXi> Cell2Ds;
          std::vector<Cell3D> Cell3Ds;
      };

      OpenVolumeMeshInterface();
      ~OpenVolumeMeshInterface();

      OVMMesh StringsToOVMMesh(const std::vector<std::string>& fileLines) const;
      std::vector<std::string> OVMMeshToStrings(const OVMMesh& mesh) const;
      OVMMesh MeshDAOToOVMMesh(const IMeshDAO& originalMesh,
                               const std::vector<std::vector<bool>>& cell3DsFacesOrientation) const;
      void OVMMeshToMeshDAO(const OVMMesh& originalMesh,
                            IMeshDAO& convertedMesh,
                            std::vector<std::vector<bool>>& convertedMeshCell3DsFacesOrientation) const;

      void ImportMeshFromFile(const std::string& ovmFilePath,
                              IMeshDAO& mesh,
                              std::vector<std::vector<bool>>& meshCell3DsFacesOrientation) const;

      void ExportMeshToFile(const IMeshDAO& mesh,
                            const std::vector<std::vector<bool>>& meshCell3DsFacesOrientation,
                            const std::string& ovmFilePath) const;

  };
}

#endif // __OpenVolumeMeshInterface_H

