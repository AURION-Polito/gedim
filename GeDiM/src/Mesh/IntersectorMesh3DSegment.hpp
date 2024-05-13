#ifndef __IntersectorMesh3DSegment_H
#define __IntersectorMesh3DSegment_H

#include "Eigen/Eigen"
#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"

namespace Gedim
{
  class IntersectorMesh3DSegment final
  {
    public:
      struct IntersectionMesh final
      {
          struct IntersectionMeshPoint final
          {
              std::vector<unsigned int> Cell2DIds = {};
              std::vector<unsigned int> Edge2DIds = {};
              std::vector<unsigned int> Vertex2DIds = {};
          };

          struct IntersectionMeshSegment final
          {
              std::vector<unsigned int> Points = {};
              std::vector<unsigned int> Cell2DIds = {};
              std::vector<unsigned int> Edge2DIds = {};
          };

          std::vector<IntersectionMeshPoint> Points;
          std::vector<IntersectionMeshSegment> Segments;
      };

    private:
      const Gedim::GeometryUtilities& _geometryUtilities;

      IntersectionMesh::IntersectionMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                                     std::map<double, IntersectionMesh::IntersectionMeshPoint>& points,
                                                                     bool& found);

    public:
      IntersectorMesh3DSegment(const Gedim::GeometryUtilities& geometryUtilities);
      ~IntersectorMesh3DSegment();

      void CreateIntersectionMesh(const Eigen::Vector3d& segmentOrigin,
                                  const Eigen::Vector3d& segmentEnd,
                                  const Eigen::Vector3d& segmentTangent,
                                  const Eigen::Vector3d& segmentBarycenter,
                                  const double& segmentLength,
                                  const Gedim::IMeshDAO& mesh3D,
                                  IntersectionMesh& result);
  };
}

#endif
