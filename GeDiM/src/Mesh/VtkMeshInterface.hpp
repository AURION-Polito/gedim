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
              enum struct Types
              {
                Unknown = 0,
                Tetrahedron = 1,
                Hexahedron = 2
              };

              std::vector<unsigned int> Cell0D_Id;
              Types Type;
          };

          unsigned int NumCell0Ds;
          unsigned int NumCell3Ds;

          Eigen::MatrixXd Cell0Ds;
          std::vector<Cell3D> Cell3Ds;
      };

      VtkMeshInterface();
      ~VtkMeshInterface();

      void VtkMeshToMeshDAO(const VtkMesh& originalMesh,
                            IMeshDAO& convertedMesh) const;

      void ImportMeshFromFile(const std::string& vtkFilePath,
                              IMeshDAO& mesh) const;

  };
}

#endif // __OpenVolumeMeshInterface_H

