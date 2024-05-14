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
      struct IntersectionPoint final
      {
          std::set<unsigned int> Cell3DIds;
      };

      struct FindSegmentStartingCell3DResult final
      {
          unsigned int StartingCell3DIndex;
          bool FoundOtherIntersections;
      };

      const Gedim::GeometryUtilities& geometryUtilities;
      const Gedim::MeshUtilities& meshUtilities;

      IntersectionPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                               std::map<double, IntersectionPoint>& points,
                                               bool& found) const;

      void CheckSegmentIntersection(const Gedim::GeometryUtilities::PointSegmentPositionTypes segment_intersection_position,
                                    const double& segment_intersection_coordinate,
                                    const Gedim::IMeshDAO& mesh3D,
                                    const unsigned int cell3D_index,
                                    const unsigned int cell2D_index,
                                    std::map<double, IntersectionPoint>& mesh1D_intersections,
                                    std::list<unsigned int>& cell3Ds_index) const;

      std::vector<IntersectionMesh::IntersectionMeshSegment> CreateIntersectionSegments(const std::vector<IntersectionMesh::IntersectionMeshPoint>& mesh1D_points) const;

      FindSegmentStartingCell3DResult FindSegmentStartingCell3D(const Eigen::Vector3d& segmentOrigin,
                                                                const Eigen::Vector3d& segmentEnd,
                                                                const Gedim::IMeshDAO& mesh3D,
                                                                const Gedim::MeshUtilities::MeshGeometricData3D& mesh3D_geometricData) const;
      std::vector<IntersectionMesh::IntersectionMeshPoint> FindSegmentIntersectionPoints(const Eigen::Vector3d& segmentOrigin,
                                                                                         const Eigen::Vector3d& segmentEnd,
                                                                                         const Eigen::Vector3d& segmentTangent,
                                                                                         const Eigen::Vector3d& segmentBarycenter,
                                                                                         const double& segmentLength,
                                                                                         const Gedim::IMeshDAO& mesh3D,
                                                                                         const Gedim::MeshUtilities::MeshGeometricData3D& mesh3D_geometricData,
                                                                                         const unsigned int starting_cell3D_index) const;

      void SegmentCell3DIntersection(const Eigen::Vector3d& segmentOrigin,
                                     const Eigen::Vector3d& segmentEnd,
                                     const Eigen::Vector3d& segmentTangent,
                                     const Eigen::Vector3d& segmentBarycenter,
                                     const double& segmentLength,
                                     const Gedim::IMeshDAO& mesh3D,
                                     const Gedim::MeshUtilities::MeshGeometricData3D& mesh3D_geometricData,
                                     const unsigned int cell3D_index,
                                     std::map<double, IntersectionPoint>& mesh1D_intersections,
                                     std::list<unsigned int>& cell3Ds_index) const;

    public:
      IntersectorMesh3DSegment(const Gedim::GeometryUtilities& geometryUtilities,
                               const Gedim::MeshUtilities& meshUtilities);
      ~IntersectorMesh3DSegment();

      static std::vector<double> ToCurvilinearCoordinates(const IntersectorMesh3DSegment::IntersectionMesh& intersectingMesh);
      static std::string ToString(const IntersectorMesh3DSegment::IntersectionMesh& intersectingMesh);

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
