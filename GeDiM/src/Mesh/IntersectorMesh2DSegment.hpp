#ifndef __IntersectorMesh2DSegment_H
#define __IntersectorMesh2DSegment_H

#include "Eigen/Eigen"
#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"

using namespace std;

namespace Gedim
{
  class IntersectorMesh2DSegment final
  {
    public:
      struct IntersectionMesh final
      {
          struct IntersectionMeshPoint final
          {
              set<unsigned int> Cell2DIds = {};
              set<unsigned int> Edge2DIds = {};
              set<unsigned int> Vertex2DIds = {};
          };

          struct IntersectionMeshSegment final
          {
              vector<double> Points = {};
              set<unsigned int> Cell2DIds = {};
              set<unsigned int> Edge2DIds = {};
          };

          map<double, IntersectionMeshPoint> Points;
          vector<IntersectionMeshSegment> Segments;
      };

      /// \brief convert IntersectionMesh to Curvilinear Coordinates vector
      static void ToCurvilinearCoordinates(const IntersectionMesh& intersectingMesh,
                                           vector<double>& curvilinearCoordinates);

      static void ToString(const IntersectionMesh& intersectingMesh);

    private:
      const Gedim::IMeshDAO& _mesh;
      const Gedim::GeometryUtilities& _geometryUtilities;

      IntersectionMesh::IntersectionMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                                     IntersectionMesh& result,
                                                                     bool& found);
      void CheckOriginAndEndSegmentPosition(const Eigen::Vector3d& segmentOrigin,
                                            const Eigen::Vector3d& segmentEnd,
                                            IntersectionMesh& result);
      void CreateIntersectionPoints(const Eigen::Vector3d& segmentOrigin,
                                    const Eigen::Vector3d& segmentEnd,
                                    const Eigen::Vector3d& segmentTangent,
                                    const Eigen::Vector3d& segmentBarycenter,
                                    const double& segmentLength,
                                    IntersectionMesh& result);
      void CreateIntersectionSegments(IntersectionMesh& result);

      void CheckOriginAndEndPointExistence(IntersectionMesh& result);

      void SmoothIntersections(const Eigen::Vector3d& segmentOrigin,
                               const Eigen::Vector3d& segmentEnd,
                               IntersectionMesh& result);

    public:
      IntersectorMesh2DSegment(const Gedim::IMeshDAO& mesh,
                               const Gedim::GeometryUtilities& geometryUtilities);
      ~IntersectorMesh2DSegment();

      void CreateIntersectionMesh(const Eigen::Vector3d& segmentOrigin,
                                  const Eigen::Vector3d& segmentEnd,
                                  const Eigen::Vector3d& segmentTangent,
                                  const Eigen::Vector3d& segmentBarycenter,
                                  const double& segmentLength,
                                  IntersectionMesh& result);
  };
}

#endif
