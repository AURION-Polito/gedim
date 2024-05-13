#include "IntersectorMesh3DSegment.hpp"
#include "IOUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  IntersectorMesh3DSegment::IntersectorMesh3DSegment(const Gedim::GeometryUtilities& geometryUtilities) :
    _geometryUtilities(geometryUtilities)
  {
  }
  IntersectorMesh3DSegment::~IntersectorMesh3DSegment()
  {
  }
  // ***************************************************************************
  IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshPoint& IntersectorMesh3DSegment::InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                                                     std::map<double, IntersectionMesh::IntersectionMeshPoint>& points,
                                                                                                                     bool& found)
  {
    double foundCoordinate = -1.0;
    for (auto& it : points)
    {
      if (!_geometryUtilities.IsValuePositive(abs(it.first - curvilinearCoordinate),
                                              _geometryUtilities.Tolerance1D()))
      {
        foundCoordinate = it.first;
        break;
      }
    }

    if (foundCoordinate != -1.0)
    {
      found = true;
      return points[foundCoordinate];
    }

    points.insert(std::make_pair(curvilinearCoordinate,
                                 IntersectionMesh::IntersectionMeshPoint()));
    found = false;
    return points[curvilinearCoordinate];
  }
  // ***************************************************************************
  void IntersectorMesh3DSegment::CreateIntersectionMesh(const Eigen::Vector3d& segmentOrigin,
                                                        const Eigen::Vector3d& segmentEnd,
                                                        const Eigen::Vector3d& segmentTangent,
                                                        const Eigen::Vector3d& segmentBarycenter,
                                                        const double& segmentLength,
                                                        const Gedim::IMeshDAO& mesh3D,
                                                        IntersectionMesh& result)
  {

  }
  // ***************************************************************************
}
