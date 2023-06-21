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
          unsigned int NumCell0Ds;
          unsigned int NumCell1Ds;
          unsigned int NumCell2Ds;
          unsigned int NumCell3Ds;

          Eigen::MatrixXd Cell0Ds;
          Eigen::MatrixXi Cell1Ds;
          std::vector<Eigen::MatrixXi> Cell2Ds;
          std::vector<std::vector<unsigned int>> Cell3Ds;
      };

      OpenVolumeMeshInterface();
      ~OpenVolumeMeshInterface();

      OVMMesh ConvertOVMMesh(const std::vector<std::string>& fileLines) const;
      void ConvertGedimMesh(const OVMMesh& meshImported,
                            IMeshDAO& convertedMesh) const;

      void ImportMesh(IMeshDAO& mesh,
                      const std::string& ovmFilePath) const;

  };
}

#endif // __OpenVolumeMeshInterface_H

