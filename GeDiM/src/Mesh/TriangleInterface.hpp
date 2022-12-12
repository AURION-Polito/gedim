#ifndef __TriangleInterface_H
#define __TriangleInterface_H

#include "Eigen/Eigen"
#include "mytriangle.hpp"
#include "IMeshDAO.hpp"

namespace Gedim
{
  class TriangleInterface final
  {
    private:
      void CreateTriangleInput(const Eigen::MatrixXd& polygonVertices,
                               struct triangulateio& triangleInput,
                               const Eigen::MatrixXd& constrainedPoints = Eigen::MatrixXd(),
                               const Eigen::MatrixXi& constrainedSegments = Eigen::MatrixXi()) const;
      void CreateTriangleOutput(const double& maxTriangleArea,
                                struct triangulateio& triangleInput,
                                struct triangulateio& triangleOutput,
                                const std::string& triangleOptions = "-QDzpqnea") const;

      void DeleteTriangleStructure(struct triangulateio& triangleInput,
                                   struct triangulateio& triangleOutput) const;

      void ExportTriangleMesh(const struct triangulateio& triangleInput,
                              const struct triangulateio& triangleOutput,
                              const std::string& nameFolder,
                              const std::string& nameFile) const;
    public:
      TriangleInterface();
      ~TriangleInterface();

      void CreateMesh(const Eigen::MatrixXd& polygonVertices,
                      const double& maxTriangleArea,
                      IMeshDAO& mesh,
                      const std::string& triangleOptions = "-QDzpqnea") const;

  };
}

#endif // __TriangleInterface_H

