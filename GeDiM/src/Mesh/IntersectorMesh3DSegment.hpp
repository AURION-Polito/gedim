#ifndef __IntersectorMesh3DSegment_H
#define __IntersectorMesh3DSegment_H

#include "Eigen/Eigen"
#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
  class IntersectorMesh3DSegment final
  {
    public:
      struct IntersectionMesh final
      {
          struct IntersectionMeshPoint final
          {
              double CurvilinearCoordinate;
              std::vector<unsigned int> Cell3DIds;
          };

          struct IntersectionMeshSegment final
          {
              std::array<unsigned int, 2> PointsIndex;
              std::vector<unsigned int> Cell3DIds;
          };

          std::vector<IntersectionMeshPoint> Points;
          std::vector<IntersectionMeshSegment> Segments;
      };

    private:
      const Gedim::GeometryUtilities& geometryUtilities;
      const Gedim::MeshUtilities& meshUtilities;

      IntersectionMesh::IntersectionMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                                     std::map<double, IntersectionMesh::IntersectionMeshPoint>& points,
                                                                     bool& found) const;

      std::vector<IntersectionMesh::IntersectionMeshSegment> CreateIntersectionSegments(const std::vector<IntersectionMesh::IntersectionMeshPoint>& mesh1D_points) const;

      unsigned int FindSegmentVertexCell3D(const Eigen::Vector3d& vertex,
                                           const Gedim::IMeshDAO& mesh3D,
                                           const Gedim::MeshUtilities::MeshGeometricData3D& mesh3D_geometricData) const;

    public:
      IntersectorMesh3DSegment(const Gedim::GeometryUtilities& geometryUtilities,
                               const Gedim::MeshUtilities& meshUtilities);
      ~IntersectorMesh3DSegment();

      IntersectionMesh CreateIntersectionMesh(const Eigen::Vector3d& segmentOrigin,
                                              const Eigen::Vector3d& segmentEnd,
                                              const Eigen::Vector3d& segmentTangent,
                                              const Eigen::Vector3d& segmentBarycenter,
                                              const double& segmentLength,
                                              const Gedim::IMeshDAO& mesh3D,
                                              const Gedim::MeshUtilities::MeshGeometricData3D& mesh3D_geometricData) const;
  };
}

#endif
