#ifndef __VtkMeshInterface_H
#define __VtkMeshInterface_H

#include "IMeshDAO.hpp"

namespace Gedim
{
  class VtkMeshInterface final
  {
    public:
      struct VtkMesh final
      {
          struct Cell3D final
          {
              std::vector<unsigned int> FacesIndex;
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

      VtkMeshInterface();
      ~VtkMeshInterface();

      VtkMesh StringsToVtkMesh(const std::vector<std::string>& fileLines) const;
      std::vector<std::string> VtkMeshToStrings(const VtkMesh& mesh) const;
      VtkMesh MeshDAOToVtkMesh(const IMeshDAO& originalMesh) const;
      void VtkMeshToMeshDAO(const VtkMesh& originalMesh,
                            IMeshDAO& convertedMesh) const;

      void ImportMeshFromFile(const std::string& vtkFilePath,
                              IMeshDAO& mesh) const;

      void ExportMeshToFile(const IMeshDAO& mesh,
                            const std::string& vtkFilePath) const;

  };
}

#endif // __OpenVolumeMeshInterface_H

