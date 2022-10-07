#ifndef __TetgenInterface_H
#define __TetgenInterface_H

#include "IMeshDAO.hpp"
#include "tetgen.h"
#include "Eigen/Eigen"

namespace Gedim
{
  /// \brief The Tetgen Interface
  /// \see http://wias-berlin.de/software/tetgen/files/tetcall.cxx
  class TetgenInterface final
  {
    protected:

      std::string tetgenOptions;


      void CreateTetgenInput(const Eigen::MatrixXd& polyhedronVertices,
                             const Eigen::MatrixXi& polyhedronEdges,
                             const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                             tetgenio& tetgenInput,
                             const Eigen::MatrixXd& constrainedPoints = Eigen::MatrixXd(),
                             const std::vector<Eigen::VectorXi>& constrainedFaces = std::vector<Eigen::VectorXi>()) const;
      void CreateTetgenOutput(const double& maxTetrahedronArea,
                              tetgenio& tetgenInput,
                              tetgenio& tetgenOutput,
                              const std::string& tetgenOptions = "Qpqfezna") const;
      void ConvertTetgenOutputToMeshDAO(const tetgenio& tetgenOutput,
                                        IMeshDAO& mesh) const;

      void DeleteTetgenStructure(tetgenio& tetgenInput,
                                 tetgenio& tetgenOutput) const;

    public:
      TetgenInterface();
      ~TetgenInterface();

      void CreateMesh(const Eigen::MatrixXd& polyhedronVertices,
                      const Eigen::MatrixXi& polyhedronEdges,
                      const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                      const double& maxTetrahedronVolume,
                      IMeshDAO& mesh) const;

      void ExportTetgenOutput(const std::string& nameFolder,
                              const std::string& nameFile,
                              tetgenio& tetgenOutput) const;
  };

}

#endif // __TetgenInterface_H
